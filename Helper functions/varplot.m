function [ handle ] = varplot( varargin )
%VARPLOT: Plot data with variance indicators, in pretty Nature figure 
%style. Supports NaN values in data (treated as missing data). Can be used 
%and customized just like the standard "plot" function - for example: 
%varplot(data,'r','linewidth',2) etc.. Can plot confidence intervals, 
%standard deviation, and standard error (default), or use custom shading
%borders. 
% 
%Usage: 
%VARPLOT(D), where D is a matrix, plots the mean and variance of the
%matrix (the mean is plotted as a line and the variance as a shading around
%it). The data is assumed to be arranged in columns, such that each column 
%is a "trial". The default variance indicator is 95% confidence intervals. 
%
%VARPLOT(x,D), where x is a vector of the same length as D, defines the 
%x-axis of the plot, as in the plot function. 
%
%VARPLOT(v,low,high), where v, low and high are vectors of equal length,
%plots v with the shading borders defined by the low and high vectors. 
%
%VARPLOT(...,method), where method is either 'ci','std' or 'stderr',
%defines the variance indicator as confidence intervals, one standard
%deviation, or one standard error from the mean. 'ci' must be followed by 
%the confidence interval alpha, for instance VARPLOT(D,'ci',0.9). Default
%method is standard error. 
%
%VARPLOT(...,'transparency',t), where t is between 0 and 1, defines the 
%transparency level of the shading. Default is 0.5. If a value of 1 is
%given, dotten lines will be plotted instead of the transparent patch.
%
%VARPLOT(...,'palefactor',p), where p is between 0 and 1, defines how
%much paler (whiter) the shading will be compared to the line (in addition
%to the transparency effect). Setting the pale factor to 0 will use the 
%same color, while setting it to 1 will always change the color to white.
%Default is 0.5.
%
%VARPLOT(...,<plot arguments>) will pass additional arguments to Matlab's 
%standard plot function. Note that these arguments must always come after 
%the varplot-specific arguments. For instance, VARPLOT(D,'transparency',
%0.2,'k','linewidth',2) will set the line color and line width in the same 
%way as it would when using the plot function. 
% 
%h = VARPLOT(...) will return a handle to the plot object, in the same way
%as for the plot function. 
%
% 
% Written by Edden M. Gerber, lab of Leon Y. Deouell, April 2014. 
% Please send bug reports and requsts to edden.gerber@gmail.com
%
% 
% CHANGE LOG:
% 
% 20/5/2015 Edden Gerber: Added support for NaN values. 1) Variance
%           shading will now be added using separate patches between NaN 
%           values. 2) CI, std.dev. and std.err. are now ignoring missing 
%           values for time points where only part of the values are NaN. 
% 20/5/2015 Edden Gerber: Changed the default method from CI 95% to
%           std.err.
% 28/7/2015 Edden Gerber: added a "uistack" call so that lines are visually
%           displayed on top of other objects. Fixed the problem this
%           creates with the legend order by directly setting the 
%           "LegendInformation" property of the patches so they would be
%           ignored by the legend. Removed the earlier fix which relied on
%           the object order. 
% 10/4/2016 Corrected a bug in the "findobj" line at the end, which
%           searched for line objects across all matlab figures instead of
%           just in the current axis (this caused run time to increase when
%           there were many graphical objects overall). 
%
% 13/7/2016 Carmel Asch: added verification that low and high limits (for 
%           vector input data) are column vectors. The data is converted to
%           column form and if limits are not it produces an error in "fill" 
%           function - section of adding shade. My addition was performed 
%           in section that defines y vector.

% Handle input
var_method = 'stderr';
ci_alpha = 0.95;
transparency = 0.5;
pale_factor = 0.5;
ci_meanVal=0;

if nargin < 1
    error('No input arguments');
end

% count number of initial numeric arguments
num_numeric = 0;
for ind = 1:nargin
    if isnumeric(varargin{ind})
        num_numeric = num_numeric + 1;
    else
        break;
    end
end

switch num_numeric
    case 0
        error('First input argument must be vector or matrix.');
    case 1
        data = varargin{1};
    case 2
        x_vec = varargin{1};
        data = varargin{2};
    case 3
        data = varargin{1};
        lo_lim = varargin{2};
        hi_lim = varargin{3};
    case 4 
        x_vec = varargin{1};
        data = varargin{2};
        lo_lim = varargin{3};
        hi_lim = varargin{4};
    otherwise
        error('Not expecting more than five initial numeric input arguments');
end

% if data is a vector, make sure it's a column vector
if size(data,1)==1
    data = data';
end

% make an x vector variable if it was not defined.
if ~exist('x_vec','var')
    x_vec = 1:size(data,1);
end
if size(x_vec,1)==1;
    x_vec = x_vec';
end

% define the y vector
if isvector(data)
    matrix_input = false;
    y_vec = data;
    if ~exist('lo_lim','var') || ~exist('hi_lim','var')
        error('Non-matrix input must be followed by lower-limit vector and upper-limit vector');
    end
    if length(lo_lim) ~= length(data) || length(hi_lim) ~= length(data)
        error('Lower-limit vector, upper-limit vector and data vector must be of equal length');
    end
    % Make lo_lim and hi_lim are column vectors like data
    if size(lo_lim,1)==1
        lo_lim = lo_lim';
    end
    if size(hi_lim,1)==1
        hi_lim = hi_lim';
    end
else
    matrix_input = true;
    y_vec = nanmean(data,2);
end

% handle additional optional input arguments
plot_arguments = {};
ind = num_numeric + 1;
while ind <= nargin
    switch varargin{ind}
        case 'ci'
            if ind == nargin | ~isnumeric(varargin{ind+1})
                error('Optional argument ''ci'' must be followed by the confidence interval alpha, for example 0.95');
            end
            var_method = 'ci';
            ci_alpha = varargin{ind+1};
            ind = ind + 2;
        case 'std'
            var_method = 'std';
            ind = ind + 1;
        case 'stderr'
            var_method = 'stderr';
            ind = ind + 1;
        case 'palefactor'
            if ~isnumeric(varargin{ind+1}) || varargin{ind+1} < 0 || varargin{ind+1} > 1
                error('Optional argument ''palefactor'' must be followed by a numeric between 0 and 1');
            end
            pale_factor = varargin{ind + 1};
            ind = ind + 2;
        case 'transparency'
            if ~isnumeric(varargin{ind+1}) || varargin{ind+1} < 0 || varargin{ind+1} > 1
                error('Optional argument ''transparency'' must be followed by a numeric between 0 and 1');
            end
            transparency = varargin{ind + 1};
            ind = ind + 2;
        case 'ci_ttestMeanValue'
            ci_meanVal = varargin{ind + 1};
            ind = ind + 2;
        otherwise
            plot_arguments = varargin(ind:end);
            break;
    end
end

% Determine shade limits
if matrix_input
    switch var_method
        case 'ci'
            [~,~,ci] = ttest(data',ci_meanVal,'alpha',1-ci_alpha);
            lo_lim = ci(1,:)';
            hi_lim = ci(2,:)';
        case 'std'
            stdev = nanstd(data,0,2);
            lo_lim = y_vec - stdev;
            hi_lim = y_vec + stdev;
        case 'stderr'
            stderr = nanstd(data,0,2) ./ sqrt(sum(~isnan(data),2));
            lo_lim = y_vec - stderr;
            hi_lim = y_vec + stderr;
    end
end

% Plot mean
plot_held = ishold;
h_plot = plot(x_vec,y_vec,plot_arguments{:});
xlim([x_vec(1) x_vec(end)]);
    
% add shading
if transparency < 1
    clr = get(h_plot,'color');
    clr = clr + pale_factor*(1-clr);
    
    % speficiy patches - if there NaNs in the time series, they act as
    % patch borders
    nans_in_y_vec = isnan(y_vec);
    start_patch = find(diff([1;nans_in_y_vec])==-1);
    end_patch = find(diff([nans_in_y_vec;1])==1);
    
    hold all
    for p = 1:length(start_patch)
        x1 = start_patch(p);
        x2 = end_patch(p);
        h = fill([x_vec(x1:x2) ; x_vec(x2:-1:x1)],[lo_lim(x1:x2) ; hi_lim(x2:-1:x1)],clr,'FaceAlpha',1-transparency,'EdgeAlpha',0,'Tag','NotInLegend');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; % this removed the patch object from the legend list
    end
else
    hold all
    h_lo = plot(x_vec,hi_lim,plot_arguments{:});
    h_hi = plot(x_vec,lo_lim,plot_arguments{:});
    set(h_lo,'LineStyle','--','Color',get(h_plot,'Color'),'LineWidth',get(h_plot,'LineWidth')/2,'Tag','NotInLegend'); % same color, half width, dotted
    set(h_hi,'LineStyle','--','Color',get(h_plot,'Color'),'LineWidth',get(h_plot,'LineWidth')/2,'Tag','NotInLegend'); % same color, half width, dotted
end
% if the plot was not help previously, keep it unheld
if ~plot_held 
    hold off
end

% made obsolete by directly setting the LegendInformation property of the
% patch objects. 
%{
% re-sort the axis children obejct handles, so that the first handles will
% refer to the lines and not the patches (this is necessary so that the 
% lines will appear in a legend and not the patches). 
child_handles = get(gca,'Children');
tags = get(child_handles,'Tag');
not_in_legend = cellfun(@strcmp,tags,repmat({'NotInLegend'},length(tags),1));
set(gca,'Children',[child_handles(not_in_legend) ; child_handles(~not_in_legend)]);
%}

% re-sort the visual stacking order such that lines are on top (not hidden
% under shading patches)
line_handles = findobj(gca,'Type','Line');
for i=length(line_handles):-1:1
    uistack(line_handles(i),'top');
end

% output plot handle
if nargout > 0
handle = h_plot;

end

