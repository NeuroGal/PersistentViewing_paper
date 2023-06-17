function fig = plot_results_mat(result_matrix, time_vec, title_str, varargin) 
% Written by Gal Vishne, Deouell Lab 2021-2022
% Bug reports \ requests: gal.vishne@gmail.com
%
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% It can be used to plot simple line plots, matrix info (imagesc) and even
%   matrix + 'slices' of it as line plots (see Figure 3d of the paper).
% The function was written in multiple batches, so sorry it is a bit messy,
%   I think you can still use it well if you want to reproduce things I did
%   in the paper - and you are free to contact me with questions.
%
% * uses:
%    - https://colorbrewer2.org/ by Cynthia Brewer and Mark Harrower
%    - linspecer: https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap
%    - Tight_subplot: https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
% ** linspecer and cbrewer were replaced by the RGB values here so you
% don't need to download them. Tight_subplot is in the folder.
%
% Inputs:
%   result_matrix: eg rdm correlations between RDMs n_time x n_time or
%               temporal generalization decoding matrix
%              - I also incorporated plotting a single line (n_time x 1) and
%               you have my function nice_line_plot for that.
%              - It can also get a matrix that's not symmetric (n_time1 x n_time 2). 
%                That was useful for me, not sure if you need it.
%   time_vec: corresponding to n_time in ms (can be cell if rows & cols have
%               different time, see below)
%   title_str: made it mandatory so everything is clearly marked (:
%
% optional:
% stats:
%   'mask'            - mark a region of the plot, e.g. because it's
%                       significant. has two modes - either mark with
%                       surrounding lines, or fade the rest (using
%                       'FaceAlpha'). currently it's fixed to 'lines'
%                       below.
%                       If plotting line input it will add a horizontal bar
%                       to mark significant points. That's also how the
%                       'time slices' are marked (see below). They are just
%                       based on the stats of the main matrix (also if it's
%                       plotting the diagonal).
%   'mask_info'       - visualization for the mask line, assumed to be followed by
%                       struct with fields  - color,offset (see below units)
%   'cluster_pval'    - followed by array with the pvalues (relevant
%                       for cluster stats, added as text). currently I made
%                       it default to say what all clusters pass, but can
%                       be switched to largest by using another input
%   'largest_clust'   - present the p value of the largest cluster

% color input:
%   'colormap'        - default: CT=cbrewer('seq', 'YlGnBu', 64, 'spline');
%                       % can be just [r g b] for vector results.
%   'line_col'        - followed by 1x3 vec (for the slices) & also used to
%                       color the titles
% 
%   'time_slices'     - slices of the big matrix to present on the right.
%       insert as: {[t1 t2], [t3 t4], ...}, then 1st plot will be average
%       of rows corresponding to tp t1-t2, 2nd average of t3-t4 etc. 
%       (default: 150 : 300, time_vec(end)-150 : time_vec(end) (both in ms))
%       if empty cell doesn't plot any time slices.
%       * it also adds the diagonal like this, unless you tell it not to
%       with the next optional argument.
%   'no_diag'         - removes the diagonal from the right, or if you
%       didn't ask for diagonals then the default is to add a dashed line 
%       on the diagonal, so this removes it.
%
% visualization misc:
%   'no_fig'          - don't open new figure
%   'clims'           - limits for the plots (= ylim for the time slices)
%                       (default is according to the mat)
%   'extend_time'     - followed by t \ [t1, t2] - extend the x&y axes until
%                       this number (if just 1 assumes same for both)
%   'chance_level'    - default: 0. used to add dashed lines to mark this
%                       in the slice plots, and to adjust the colormap
%                       around it.(!)
%   'metric'          - for labels - default: 'Corr'. (should have really
%                       just called it 'ylabel')
%   'stim_len'        - in *ms* (to plot stim end line)
%   'decoding'        - changes the x\y labels - patchy, I know
%   'binning'         - followed by bin size - changes x\y lims - NVM
%
% Written by Gal Vishne, Deouell Lab 2021-2022
% Bug reports \ requests: gal.vishne@gmail.com

if ~iscell(time_vec); time_vec = repmat({time_vec},1,2); end
time_slices = {[150 300],[time_vec{1}(end)-150 time_vec{1}(end)]};
if time_vec{1}(end) < 350; time_slices = time_slices(1); end
plot_diag = true;  
if diff(cellfun(@length,time_vec))~=0; plot_diag = false; end
new_fig = true; mark_data_end = false;

mask = ones(size(result_matrix)); show_mask = false;
% CT = flipud(min(max(cbrewer('div','RdBu', 25, 'spline'),0),1)); 
CT = [0.0196078431372549,0.188235294117647,0.380392156862745;0.0980392156862745,0.317647058823529,0.572549019607843;
    0.129411764705882,0.4,0.674509803921569;0.145098039215686,0.454901960784314,0.721568627450980;
    0.176470588235294,0.505882352941176,0.741176470588235;0.262745098039216,0.576470588235294,0.764705882352941;
    0.411764705882353,0.678431372549020,0.815686274509804;0.572549019607843,0.772549019607843,0.870588235294118;
    0.682352941176471,0.831372549019608,0.901960784313726;0.752941176470588,0.862745098039216,0.917647058823529;
    0.819607843137255,0.898039215686275,0.941176470588235;0.901960784313726,0.945098039215686,0.972549019607843;
    0.968627450980392,0.968627450980392,0.968627450980392;0.992156862745098,0.929411764705882,0.890196078431373;
    0.992156862745098,0.858823529411765,0.780392156862745;0.992156862745098,0.8,0.686274509803922;
    0.984313725490196,0.737254901960784,0.6;0.956862745098039,0.647058823529412,0.509803921568627;
    0.898039215686275,0.517647058823530,0.403921568627451;0.839215686274510,0.376470588235294,0.301960784313725;
    0.796078431372549,0.258823529411765,0.235294117647059;0.760784313725490,0.164705882352941,0.192156862745098;
    0.698039215686275,0.0941176470588235,0.168627450980392;0.588235294117647,0.0392156862745098,0.149019607843137;
    0.403921568627451,0,0.121568627450980];
chance_level = 0; metric = 'Corr'; bin_size = 0;
stim_len = 0; % just used for plotting lines
mask_mode = 'lines'; % alpha\lines % if changing to alpha change the default mask to show everything
largest_clust = false;
no_colormap_given = true; decoding = false; pval_text = false;
mask_info = [];
mask_info.offset = 0.075; % from the bottom, proportion of y range
mask_info.color = 'r'; line_col = 'k';
arg  = 1;
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'time_slices'
            time_slices = varargin{arg+1};
            arg = arg + 2;
        case 'mask'
            mask = varargin{arg+1};
            show_mask = true;
            arg = arg + 2;
        case 'colormap'
            CT = varargin{arg+1};
            arg = arg + 2;
            no_colormap_given = false;
        case 'no_fig'
            new_fig = false;
            arg = arg + 1;
        case 'line_col'
            line_col = varargin{arg+1};
            arg = arg + 2;
        case 'clims'
            clims = varargin{arg+1};
            arg = arg + 2;
        case 'metric'
            metric = varargin{arg+1};
            arg = arg + 2;
        case 'chance_level'
            chance_level = varargin{arg+1};
            arg = arg + 2;
        case 'stim_len'
            stim_len = varargin{arg+1};
            if isempty(stim_len); stim_len = 0; end
            arg = arg + 2;
        case 'extend_time'
            axes_end = varargin{arg+1}; mark_data_end = true;
            if isscalar(axes_end); axes_end = repmat(axes_end, 1, 2); end
            arg = arg + 2;
        case 'mask_mode'
            mask_mode = varargin{arg+1};
            arg = arg + 2;
        case 'decoding'
            decoding = true;
            arg = arg + 1;
        case 'cluster_pval'
            pval_text = true;
            pval = varargin{arg + 1};
            arg = arg + 2;
        case 'binning'
            bin_size = varargin{arg+1};
            arg = arg + 2;
        case 'no_diag'
            plot_diag = false;
            arg = arg + 1;
        case 'mask_info'
            mask_info_input = varargin{arg+1};
            if isfield(mask_info_input, 'offset'); mask_info.offset = mask_info_input.offset; end
            if isfield(mask_info_input, 'color'); mask_info.color = mask_info_input.color; end
            arg = arg + 2;
        case 'largest_clust'
            largest_clust = true;
            arg = arg + 1;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end
if largest_clust; clust_type = 'max'; else; clust_type = 'all'; end

xlims = [time_vec{2}(1)-bin_size/2 time_vec{2}(end)+bin_size/2];
mat_ylims = [time_vec{1}(1)-bin_size/2 time_vec{1}(end)+bin_size/2];
if ~mark_data_end; axes_end = [xlims(2) mat_ylims(2)]; end
if isscalar(stim_len); stim_len = repmat(stim_len, 1, 2); end
if isvector(result_matrix) || ~new_fig; time_slices = {}; end
n_time_slices = length(time_slices);

if new_fig
    fig_sz = [0.15 0.15 0.42 0.4];
    pos = [0.14 0.13 0.41 0.79];
    if n_time_slices == 0; fig_sz(3) = 0.3; pos(3) = 0.8; end
    if isvector(result_matrix); pos([1,2,4]) = [0.1 0.2 0.6]; end
    fig = figure('Units','Normalized','Position',fig_sz);
    main_ax = axes('Units','normalized','Position',pos);
end

if ~exist('clims', 'var')
    clims = limer(results_matrix);
end
if no_colormap_given % then we can mess with it
    [CT, clims] = adjust_CT(CT, clims, chance_level, "cropped"); % really here for historic reasons but you don't need it, the function is in the folder
end
fig_col = get(gcf,'Color');

lims = struct('x',[xlims(1) axes_end(2)],'y', clims);
added_lines = struct('stim', stim_len(2), 'chance', chance_level);
make_nice = struct('title_str',title_str,'tit_change',false,'ylb',metric,'keep_ticks',true);
mask_info.show = show_mask; mask_info.mask = mask;
color_info = struct('ax_col',fig_col);

if ~isvector(result_matrix)
    imagesc_inputs = {time_vec{2}, time_vec{1}, result_matrix};
    if strcmp(mask_mode,'alpha')
        imagesc_inputs = [imagesc_inputs,{'AlphaData',mask}];
    end
    imagesc(imagesc_inputs{:});hold on
    mask(isnan(mask)) = false;
    if show_mask && strcmp(mask_mode,'lines')
        b = bwboundaries(mask);
        for k = 1:length(b)
           boundary = b{k};
           plot(time_vec{2}(boundary(:,2)), time_vec{1}(boundary(:,1)), 'k', 'LineWidth', 1.15);
        end
    end
    plot([0 0],mat_ylims,'k:'); plot(xlims,[0 0],'k:') % zero lines
    plot(xlims,[stim_len(1) stim_len(1)],'k:'); plot([stim_len(2) stim_len(2)], mat_ylims, 'k:') % stim end lines
    if mark_data_end && strcmp(mask_mode,'alpha')
        plot([xlims(2) xlims(2)],mat_ylims,'k-'); plot(xlims,[mat_ylims(2) mat_ylims(2)],'k-');
    end
    if plot_diag && n_time_slices == 0
        plot(xlims,mat_ylims,'k:')
    end
    if pval_text
        add_pval(pval,[axes_end(2) mat_ylims(1)],clust_type);
    end
    colormap(gca,CT); cb = colorbar;
    if new_fig
        set(cb,'Units','Normalized','Position',[0.04 0.3 0.02 0.4]); title(cb,metric);
    end
    if decoding
        xlabel('Test time (ms)');ylabel('Train time (ms)');
    else
        xlabel('Time (ms)'); ylabel('Time (ms)');
    end
    xlim([xlims(1) axes_end(2)]); ylim([mat_ylims(1) axes_end(1)]); caxis(clims)
    set(gca,'YDir','normal','FontSize',12.5,'Tag','matrix','Color',color_info.ax_col,'XTick',[0:300:axes_end(2)],'YTick',[0:300:axes_end(1)]); box off;
    title(title_str,'FontSize',14,'Color',line_col)
else
    if no_colormap_given
        color_info.line_col = [0.2157    0.4941    0.7216]; % linspecer(1);
    else
        color_info.line_col = CT(1,:);
    end 
    slice_plot(result_matrix, time_vec{2}, mask_info, lims, added_lines, color_info, make_nice)
    if mark_data_end && strcmp(mask_mode,'alpha')
        plot([xlims(2) xlims(2)],clims,'k-');
    end
    if pval_text
        add_pval(pval,[axes_end(2) lims.y(2)-0.15*range(lims.y)],clust_type);
    end
    title(title_str,'FontSize',14,'Color',line_col)
end

if n_time_slices > 0
if plot_diag; n_time_slices = n_time_slices + 1;end
if strcmp(metric,'acc') || strcmp(metric,'auc') || strcmp(metric,'AUC')
    lims.y(1) = max(lims.y(1),floor(min(result_matrix(:))*10)/10);
end
make_nice.tit_change = true; 
ha = tight_subplot(n_time_slices, 1, 0.05, [0.21 0.08], [0.62, 0.05]);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]
for sb = 1:n_time_slices
    axes(ha(sb))
    col = line_col;
    ts_idx = sb;
    if plot_diag && sb == 1
        n_time = length(time_vec{1});
        diag_width = 9; if n_time <=4; diag_width = 1; end
        mid_diag = ceil(diag_width/2);
        dat = nan(diag_width,n_time); dat_mask = false(diag_width, n_time);
        for i = 1:diag_width
            tmp = abs(i - mid_diag);dat_lims = (1+ceil(tmp/2)) : (n_time-floor(tmp/2));
            dat(i, dat_lims) = diag(result_matrix,i-mid_diag);
            dat_mask(i, dat_lims) = diag(mask,i-mid_diag);
        end
        dat = nanmean(dat,1); dat_mask = mean(dat_mask,1)>=0.5;
        tit_str = 'Diagonal';
        box_lims = [1,mid_diag,n_time,n_time,n_time-mid_diag+1,1;1,1,n_time-mid_diag+1,n_time,n_time,mid_diag]; box_lims = time_vec{1}(box_lims);
    else
        if plot_diag; ts_idx = ts_idx-1; end
        dat = mean(result_matrix(time_vec{1} >= time_slices{ts_idx}(1) & time_vec{1} <= time_slices{ts_idx}(2),:),1);
        dat_mask = mean(mask(time_vec{1} >= time_slices{ts_idx}(1) & time_vec{1} <= time_slices{ts_idx}(2),:),1)>=0.5;
        tit_str = sprintf('%d-%dms',time_slices{ts_idx}(1),time_slices{ts_idx}(2));
        box_lims = [reshape(repmat(xlims+0.01*range(xlims)*[-1 1],2,1),1,[]);time_slices{ts_idx},fliplr(time_slices{ts_idx})];
        ypos1 = pos(2)+pos(4)*(mean(time_slices{ts_idx})-mat_ylims(1))/(axes_end(1)-mat_ylims(1));        
        ypos2 = ha(sb).Position(2)+ha(sb).Position(4)/2;
        annotation('arrow',[pos(1)+pos(3) 0.62],[ypos1 ypos2],'LineWidth',0.35,'HeadStyle','plain','HeadLength',5,'HeadWidth',6)
    end
    box_lims = [box_lims, box_lims(:,1)];
    if sb == n_time_slices; keep_ticks = true; else; keep_ticks = false; end
    color_info.line_col = col;
    mask_info.mask = dat_mask;%false(size(dat_mask)); %
    make_nice.keep_ticks = keep_ticks; make_nice.title_str = tit_str;
    slice_plot(dat, time_vec{2}, mask_info, lims, added_lines, color_info, make_nice);
    axes(main_ax); plot(box_lims(1,:), box_lims(2,:), '-','Color',col,'LineWidth',1.3);
end
set(main_ax,'Clipping','off')
end
end

function slice_plot(dat, time_vec, mask_info, lims, added_lines, color_info, make_nice)
plot(lims.x, [added_lines.chance added_lines.chance], 'k:'); hold on
plot([0 0], lims.y, 'k:')
plot([added_lines.stim added_lines.stim], lims.y, 'k:')
plot(time_vec, dat, 'LineWidth', 2, 'Color', color_info.line_col);
if mask_info.show
    mask_info.mask(isnan(mask_info.mask)) = false;
    mask_dat = nan(size(time_vec)); mask_dat(boolean(mask_info.mask)) = min(lims.y) + mask_info.offset*range(lims.y);
    plot(time_vec, mask_dat, '-', 'LineWidth', 2, 'Color', mask_info.color);
end
xlim(lims.x);ylim(lims.y)
if make_nice.keep_ticks
    ylabel(make_nice.ylb); xlabel('Time (ms)'); set(gca,'XTick',[0:300:lims.x(2)])
else
    set(gca,'YTick',[],'XTick',[]);
end
set(gca,'FontSize',12,'Color',color_info.ax_col); box off
tit = title(make_nice.title_str); set(tit,'Units','Normalized');
if make_nice.tit_change
    set(tit,'FontSize',12,'Position',[0.005,1.005,0],'FontWeight','Normal','HorizontalAlignment','left');
end
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