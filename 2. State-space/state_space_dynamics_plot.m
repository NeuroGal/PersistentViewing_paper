function fig = state_space_dynamics_plot(dat, time_vec, ROIs, settings)
%   Function tailored specifically to produce the dynamics plots for the 
%   state space analyses (Figure 2 and Supplementary Figures 2-3)
% 
% Inputs
%   dat: cat\dur-pair x time x reg (if ndim=2 (time x reg) it's adjusted, also for the masks)
%   settings:
%       mask-related: use_mask (T\F), masks (cat\dur-pair x time x reg)
%       cat_idx:      (optional) which categories\duration pairs to plot
%       fig_name
%       appearance:  
%           x-axis:     x_axis, x_ticks (both determine what time-points are plotted and on which axes)
%           y-axis:     ylab, shared_y (related to ylimit)
%           title:      ROInames_full (the actual title), color_title (T\F, based on region color, if [] no title)
%           colors:     cat_cmap (line colors), reg_cmap (used for the title color & sometimes also the fill color)
%           peak marks: add_peak_line, peak_times (cat\dur x reg)
%           fill:       add_fill, fill_y (e.g. for CI limits\ max-statistic threshold)
  
if ~isfield(settings,'use_mask');      settings.use_mask = false; end
if ismatrix(dat); dat = permute(dat,[3 1 2]); end
if settings.use_mask && ismatrix(settings.masks); settings.masks = permute(settings.masks, [3 1 2]); end
if ~isfield(settings,'cat_idx') || isempty(settings.cat_idx); settings.cat_idx = 1:size(dat,1); end
if ~isfield(settings,'color_title');   settings.color_title = false; end
if ~isfield(settings,'add_peak_line'); settings.add_peak_line = false; end
if ~isfield(settings,'add_fill');      settings.add_fill = false; end
dat = dat(settings.cat_idx,:,ROIs);
n_reg = length(ROIs); n_cat = length(settings.cat_idx);

line_cols = settings.cat_cmap; 
if settings.cat_idx == 1; line_cols = [1 1 1]*0.1; end % this is the category selectivity case
fill_inputs = {'FaceAlpha',0.4,'EdgeColor','none'};
lims = struct('x',[time_vec(1) time_vec(end)],'y', []);
mask_info = struct('color',[],'offset',0.075);

nrow = 1; if n_reg > 4; nrow = 2; end; ncol = ceil(n_reg/nrow);
fig = figure('Units','Normalized','Position',[0.1 0.25 0.15+(0.64/4)*ncol 0.08+0.13*nrow], 'Name', settings.fig_name);
ha = tight_subplot(nrow, ncol, [0.105+0.08/nrow 0.01+0.08/ncol], [0.04+0.22/nrow 0.015+0.12/nrow], [0.07+0.54/ncol 0.006+0.036/ncol]);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]

all_ylims = get_ylims(squeeze(num2cell(dat,[1 2])), settings.shared_y);
for reg_i = ROIs
    axes(ha(ROIs==reg_i));
    lims.y = all_ylims(ROIs == reg_i,:);
    if settings.offset_lines
        for dur = 1:n_cat
            plot([settings.offset_times(dur) settings.offset_times(dur)],lims.y,'Color',line_cols(settings.cat_idx(dur),:),'LineWidth',0.8); hold on
        end
    end
    for cat_i = settings.cat_idx
        if n_cat == 1
            fill_col = settings.reg_cmap(reg_i,:); mask_info.color = 'r';
        else
            fill_col = line_cols(cat_i,:); mask_info.color = line_cols(cat_i, :);
        end
        if settings.add_fill
            fill_x = [time_vec,fliplr(time_vec)];
            fill_y = settings.fill_y{cat_i, reg_i};
            fill_x = fill_x(~isnan(fill_y));fill_y = fill_y(~isnan(fill_y));
            fill(fill_x,fill_y,fill_col,fill_inputs{:}); hold on;
        end
        plot(lims.x, [0 0], 'k:'); hold on; plot([0 0], lims.y, 'k:')
        plot(time_vec, squeeze(dat(settings.cat_idx==cat_i,:,ROIs==reg_i)), 'LineWidth', 2, 'Color', line_cols(cat_i,:));
        if settings.use_mask
            mask = settings.masks(cat_i,:, reg_i);
            mask(isnan(mask)) = false; 
            mask_dat = nan(size(time_vec));
            mask_dat(boolean(mask)) = min(lims.y) + mask_info.offset*range(lims.y);
            plot(time_vec, mask_dat, '-', 'LineWidth', 2, 'Color', mask_info.color);
        end
        title_str = settings.ROInames_full(reg_i); if isempty(settings.color_title); title_str = ''; end
        title_col = 'k'; if settings.color_title; title_col = settings.reg_cmap(reg_i,:); end
        title(title_str,'FontSize',14,'Color',title_col)
        if settings.add_peak_line
            peak_t = settings.peak_times(cat_i, reg_i); % assumes it's full (regions\categs not narrowed to what is plotted)
            plot([peak_t peak_t], [lims.y(1) -lims.y(1)], 'LineWidth', 1, 'Color', line_cols(cat_i,:));
        end
    end
    xlabel('Time (ms)'); xlim(lims.x); ylim(lims.y); box off
    if reg_i == ROIs(ncol*(nrow-1) + 1)
        ylb = ylabel(settings.ylab,'Rotation',0,'VerticalAlignment','middle');
        ylb.Position(1) = -0.725*range(time_vec);
    end
end
set(ha,'XTick',settings.x_ticks,'Color','w','FontSize',14); 
rm_ticks_idx = setdiff(1:length(ha),ncol*(nrow-1) + 1);
if settings.shared_y; set(ha(rm_ticks_idx), 'YTick',[]); end
if strcmp(settings.x_axis,'none')
    set(ha, 'XTick',[],'XLabel',[]);
elseif strcmp(settings.x_axis,'part')
    idx = (ncol*(nrow-1)+2):length(ha); if length(settings.x_ticks)>=5; idx = [ncol*(nrow-1)+1,idx]; end
    set(ha(idx),'XTick',settings.x_ticks)
    set(ha(rm_ticks_idx),'XLabel',[])
    set(ha(1:ncol*(nrow-1)),'XTick',[])
end
delete(ha((n_reg+1):end))
end

function all_ylims = get_ylims(array_per_reg, shared_y)
% the function is used to fix the y limits so that the 0 line will be matched across plots
% array_per_reg - cell n_regs x 1 
if isrow(array_per_reg); array_per_reg = array_per_reg'; end
all_vals = [array_per_reg{:}]; all_vals = all_vals(~isnan(all_vals));
if shared_y
    all_ylims = repmat(limer(all_vals), size(array_per_reg));
else
    all_ylims = real(cell2mat(cellfun(@(x) limer(x(~isnan(x))), array_per_reg,'UniformOutput',false)));
    z_ratio = max(0.85, mode(all_ylims(:,2)./range(all_ylims,2)));
    if z_ratio == 1
        all_ylims = repmat(limer(all_vals), size(array_per_reg));
    else
    corr_factor = z_ratio/(z_ratio-1);
    for l_idx = 1:size(array_per_reg,1)
        cur_z_ratio = all_ylims(l_idx,2)/range(all_ylims(l_idx,:));
        if cur_z_ratio > z_ratio
            all_ylims(l_idx,1) = all_ylims(l_idx,2)/corr_factor;
        else
            all_ylims(l_idx,2) = all_ylims(l_idx,1)*corr_factor;
        end
    end
    end
end
end