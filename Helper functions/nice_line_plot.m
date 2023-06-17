function nice_line_plot(data, time_vec, title_str, varargin) 
% * assumed to be one panel in an existing figure (not creating a new figure)
%
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
%
% Input:
%   data: n_time x 1
%   time_vec: corresponding to n_time in ms
%   title_str
%
% optional:
%   'mask'            - significant points to mark (done with horizontal
%                           bars, the name is borrowed from the 2d case
%                           where people sometimes only plot the
%                           significant region.
%   'mask_info'       - visualization for the mask line, assumed to be followed by
%                       struct with fields  - color,offset (see below for
%                       units)
%   'color'           - default: 'k' (can change if you want) * also colors the title
%   'ylims'           - set y axis limits
%   'ylab'            - set y axis label (x axis is set to 'Time (ms)')
%   'chance'          - default: 0. adds a dotted horizontal line to mark this.
%   'cluster_pval'    - followed by array with the pvalues (relevant
%                       for cluster stats, added as text). currently I made
%                       it default to say what all clusters pass, but can
%                       be switched to largest by using another input
%   'largest_clust'   - present the p value of the largest cluster
%   'extend_time'     - followed by t (in ms) - extend the x-axis to this limit.
%   'stim_len'        - in *ms* (add stim end line)
%
% Written by Gal Vishne, Deouell Lab 2022~
% Bug reports \ requests: gal.vishne@gmail.com

mark_data_end = false; show_mask = false;
color = 'k'; chance = 0; ylab = '';
pval_text = false; largest_clust = false;
mask_info = [];
mask_info.offset = 0.075; % from the bottom, proportion of y range
mask_info.color = 'r';
ylims = limer(data(:)); axes_end = time_vec(end);
arg  = 1;
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'mask'
            mask = varargin{arg+1};
            mask(isnan(mask)) = false; mask_dat = nan(size(data));
            show_mask = true;
            arg = arg + 2;
        case 'mask_info'
            mask_info_input = varargin{arg+1};
            if isfield(mask_info_input, 'offset'); mask_info.offset = mask_info_input.offset; end
            if isfield(mask_info_input, 'color'); mask_info.color = mask_info_input.color; end
            arg = arg + 2;
        case 'color'
            color = varargin{arg+1};
            arg = arg + 2;
        case 'ylims'
            ylims = varargin{arg+1};
            arg = arg + 2;
        case 'ylab'
            ylab = varargin{arg+1};
            arg = arg + 2;
        case 'chance'
            chance = varargin{arg+1};
            arg = arg + 2;
        case 'cluster_pval'
            pval_text = true;
            pval = varargin{arg + 1};
            arg = arg + 2;
        case 'largest_clust'
            largest_clust = true;
            arg = arg + 1;
        case 'extend_time'
            axes_end = varargin{arg+1};
            mark_data_end = true;
            arg = arg + 2;
        case 'stim_len'
            stim_len = varargin{arg+1};
            if isempty(stim_len); stim_len = 0; end
            arg = arg + 2;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end
if mark_data_end && ~exist('stim_len','var'); stim_len = 0; end
if largest_clust; clust_type = 'max'; else; clust_type = 'all'; end
if show_mask; mask_dat(boolean(mask)) = min(ylims) + mask_info.offset*range(ylims); end
xlims = [time_vec(1) axes_end];

hold on; plot(xlims, [chance chance], 'k:'); hold on; plot([0 0], ylims, 'k:')
if mark_data_end; plot([stim_len stim_len], ylims, 'k:'); end
plot(time_vec, data, 'LineWidth', 2, 'Color', color);
if show_mask; plot(time_vec, mask_dat, '-', 'LineWidth', 2, 'Color', mask_info.color); end
xlim(xlims);ylim(ylims)
ylabel(ylab); xlabel('Time (ms)');
set(gca,'FontSize',12,'XTick',0:300:xlims(2),'XTickLabelMode','auto','YTickLabelMode','auto'); box off
if pval_text; add_pval(pval,[xlims(2) ylims(2)-0.15*range(ylims)],clust_type); end
title(title_str,'FontSize',14,'Color',color)
end

function add_pval(pval,loc,type)
if ~exist('type','var')
    type = 'all';
end
if ~isempty(pval)
    clust_str = ''; start_str = '';
    if numel(pval)>1
        if strcmp(type,'max')
            value = min(pval); clust_str = 'largest-';
        elseif strcmp(type,'all')
            value = max(pval); start_str = 'all '; 
            if numel(pval)==2; start_str = 'both '; end
        end
    else
        value = pval;
    end
    str = sprintf('%sp_{%scluster} < %0.2g', start_str, clust_str, value);
    text(loc(1),loc(2),str,'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',8.5); 
end
end