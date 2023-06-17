%% Figures for results section 2 (Figure 2 & Supp Figure 2-3)
% also includes dimensionality reduction

clear all;clc;close all
%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));

load([DATA_FOLDER,'figure_settings.mat'])
fig_settings = [];
fig_settings.reg_cmap = reg_cmap;
fig_settings.ROInames_full = ROInames_full;
fig_settings.ROInames = ROInames;

%% Only long stimuli >= 900ms: Figure 2a, Supp Figure 2, Supp Figure 3e
% start by loading data to plot
load([DATA_FOLDER,'state_space_calc.mat']);
n_cat = length(cat_names);
fig_settings.cat_cmap = cat_cmap;

%% Figure 2a, Supp Figure 2a - multivariate trajectory plots
% this uses a toolbox I wrote over multiple years with A TON of other options,
% I leave it to you to explore it (state_space_plot boolbox).
% If you want to get something pretty similar to the trajectories in the
% paper use the settings below (note averaging within each category was done in advance).
% If you're interested in doing something similar with your data and need
% help with the functions please let me know! I'll be happy to assist.

views = [10 -20; 110 20; -55 75; -55 75; -25 80; -25 80];
plot_ids = 5:5:length(time_vec);
dim_red_inputs = {'link_plots',true,'n_dim', 3, 'pca_start_idx', find(time_vec>0,1)};

traj_settings = []; % in the current implementation this is needed to turn to trajectories mode, I might add more things implemented here in the future
viz_settings = []; viz_settings.link_plots = true; 
cat_inputs = []; cat_inputs.cat_inds = 1:n_cat; cat_inputs.cat_names = cat_names;
plot_settings = []; plot_settings.single = true; plot_settings.use_sem  = true; plot_settings.mean = false; plot_settings.std = false;
plot_patt_inputs  = {'trajectories', traj_settings, 'plot_settings', plot_settings, 'viz_settings', viz_settings, 'categories', cat_inputs};
for reg_i = 1:n_reg
[dim_red_patterns, pca_info] = do_dim_red(mean_segs(:, ROI_idx==reg_i, :), 'pca', plot_ids, dim_red_inputs{:});
fig = plot_mv_patterns(dim_red_patterns, time_vec(plot_ids), 'pca_info', pca_info, plot_patt_inputs{:});
set(fig,'Name',sprintf('Multivariate trajectories %s',ROInames(reg_i)))
title(findobj(fig.Children,'Tag','main_ax'),ROInames_full(reg_i),'Visible','on','FontSize',16,'Color',reg_cmap(reg_i,:)) 
view(views(reg_i,1),views(reg_i,2))
end

%% Supp Figure 2b,d plot speed diff\attenuation with ci

for opt = 1:2 % 1: attenuation (Figure S2b), 2: speed (Figure S2d)
fig_settings.ROIs       = 1:n_reg;
fig_settings.ci_pcnt    = 95;
if opt == 1
    fig_settings.ylim = [30 90];
    fig_settings.ylab = '% attenuation from peak';
    fig_settings.title = 'Multivariate response attenuation';
    fig_settings.zero_ln = false;
    fig_part = 'resp atten';
    dat = percent_atten(:,:,:,1); 
elseif opt == 2
    fig_settings.ylim = [-20 140];
    fig_settings.ylab = 'speed diff. (% from peak)';
    fig_settings.title = ['Difference in state transition speed' newline '(before peak - after peak)'];
    fig_settings.zero_ln = true;
    fig_part = 'speed diff';
    dat = speed_diff;
end
fig_settings.fig_name = sprintf('Multivariate %s ci-cat-reg', fig_part);
fig = ci_cat_reg_plot(dat, fig_settings, @get_ci_n_val);
end

%% calculation interlude: Extract range of peak\attenuation values across categories (used for Supp Figure S2 text)

meas = 1; % 1: peak_time, 2: percent_atten
dim  = 2; % 1: baseline_dist, 2: speed
if meas == 1; res_array = peak_time; elseif meas == 2; res_array = percent_atten; end
actual_val = res_array(:,:,1,dim);

minmax_vals = round([min(actual_val);max(actual_val)],1);
minmax_tbl = array2table(minmax_vals,'RowNames',{'Min','Max'},'VariableNames',ROInames(1:n_reg))
join(string(minmax_vals(1:2,:))','-')'+"ms"

%% Supp Figure 2c plot dynamics with ci (can also use opt == 1 to plot the insets for Figure 2a and Figure S2a)

opt  = 2; % 1: baseline_dist (part of Figure 2a,Figure S2a - also comes out of that function), 2: speed (Figure S2c)
ROIs = 1:n_reg;
if opt == 1; dat = baseline_dist; elseif opt == 2; dat = speed; end
dat = dat - nanmean(dat(:,time_vec<0,:,:),2); % baseline just for the visualization

% some visualization settings
fig_settings.x_axis        = 'part'; % none\part\full
fig_settings.x_ticks       = [0 900]; 
fig_settings.color_title   = true; % use [] for no title
fig_settings.use_mask      = false;

fig_settings.add_fill   = true; ci_pcnt = 95;
if fig_settings.add_fill
[~, low_ci, high_ci] = get_ci_n_val(dat(:,:,ROIs,:), 4, ci_pcnt);
fig_settings.fill_y = cellfun(@(low,high) [low, fliplr(high)], num2cell(low_ci, 2), num2cell(high_ci, 2),'UniformOutput', false);
end
fig_settings = ci_dynamics_settings(opt, false, fig_settings, ROIs, peak_time); %2nd arg = not mult_dur_mode
fig = state_space_dynamics_plot(dat(:,:,:,1), time_vec, ROIs, fig_settings);

%% Supp Figure 3e plot multivariate category selectivity index

ROIs     = 1:n_reg;
fig_settings.masks = catsel_masks;
signif_times = show_signif(fig_settings.masks, time_vec, ROInames_full); % useful for the legend etc. (tells us what are the ranges where it's significant)

fill_y = cell(1,length(ROIs));
for reg_i = ROIs
    low_ci = -100*ones(size(time_vec')); % assumes we took a one-sided test so there is no lower bound
    high_ci = catsel_thresh(reg_i)*ones(size(time_vec'));
    fill_y{ROIs==reg_i} = [low_ci; flipud(high_ci)]';
end
fig_settings.fill_y   = fill_y;
fig_settings.use_mask = true;
fig_settings.add_fill = true;

% some visualization settings
fig_settings.shared_y       = false;
fig_settings.x_axis         = 'part'; % none\part\full
fig_settings.x_ticks        = [0 900];
fig_settings.color_title    = true; % use [] for no title
fig_settings.offset_lines   = false;
fig_settings.add_peak_line  = false;
fig_settings.fig_name       = sprintf('Multivar cat var shuff-dynamics %s', join(ROInames(ROIs)));
fig_settings.ylab           = ['Multivariate category' newline 'selectivity index (z-score)'];
fig = state_space_dynamics_plot(multivar_catsel, time_vec, ROIs, fig_settings);

%% Multiple durations: Figure 2b, Supp Figure 3a-d
% start by loading data to plot
load([DATA_FOLDER,'state_space_calc_mult_dur.mat']);
n_cat = length(cat_names); n_dur = length(stim_durs); n_pair = size(dur_pairs,1);
fig_settings.cat_cmap = cat_dur_cmap;

%% plot

clc
opt         = 1; % 1 - baseline_dist (Figure 2b, Figure S3a-b),3 - dist_bet_durs (Figure S3c-d) [2 - speed but not used here]
cat_plot    = 1; % 1 - faces (Figure 2b, Figure S3a,c), 2 - watches (Figure S3b,d)
ROIs        = 1:6; % Figure 2b: 1:4, Figure S3a: 5:6, Figure S3b-d: 1:6

if opt == 1 
    dat = baseline_dist; fig_settings.masks = baseline_dist_diff_masks; % masks in both cases from the shuffled stats !
    dat = dat - mean(dat(:,time_vec<0,:),2); % baseline just for visualization
    fig_settings.cat_idx = (1:n_dur) + (cat_plot-1)*n_dur;
elseif opt == 3
    dat = dist_bet_durs; fig_settings.masks = dist_bet_durs_masks;
    fig_settings.cat_idx = (1:n_pair) + (cat_plot-1)*n_dur; % specific for the duration pairs I defined above
end
% useful for the legend etc. (tells us what are the ranges where it's significant)
signif_times = show_signif(fig_settings.masks(fig_settings.cat_idx,:,ROIs), time_vec, ROInames_full(ROIs), join(string(pair_times),'-vs-'));

% some visualization settings
fig_settings.x_axis        = 'part'; % none\part\full
fig_settings.x_ticks       = [0 300:600:2100]; 
fig_settings.color_title   = true;
fig_settings.use_mask      = true;
fig_settings.add_fill      = false;

fig_settings = ci_dynamics_settings(opt, true, fig_settings, ROIs, [], stim_durs); %2nd arg = yes mult_dur_mode
fig = state_space_dynamics_plot(dat, time_vec, ROIs, fig_settings);

%% functions

function [actual_val, low_ci, high_ci] = get_ci_n_val(array, perm_dim, ci_pcnt)
% eliminates the perm dim
n_dim = ndims(array);
array = permute(array, [1:(perm_dim-1), (perm_dim+1):n_dim, perm_dim]); % move perm dim to last
dim_sizes = size(array); n_elem = prod(dim_sizes((1:end-1)));
actual_val = reshape(array(1:n_elem), dim_sizes(1:end-1));
perm_vals = reshape(array((n_elem+1):end), [dim_sizes(1:end-1),dim_sizes(end)-1]);
low_ci  = prctile(perm_vals, (100-ci_pcnt)/2, n_dim);
high_ci = prctile(perm_vals, 100-(100-ci_pcnt)/2, n_dim);
end

function settings = ci_dynamics_settings(opt, is_mult_dur, settings, ROIs, peak_time, stim_durs)
% opt: 1 - baseline_dist, 2 - speed, 3 - dist bet durs
settings.shared_y      = false;
if exist('peak_time','var') && ~isempty(peak_time)
    settings.add_peak_line  = true;
    settings.peak_times     = peak_time(:,:,1,opt);
else
    settings.add_peak_line  = false;
end
if is_mult_dur
    type = 'duration';
    settings.offset_lines  = true;
    settings.offset_times  = stim_durs;
else
    type = 'category';
    settings.offset_lines  = false;
end
if opt == 1
    settings.ylab   = ['Global multivariate response' newline '(per ' type ', a.u.)'];
    fig_start       = 'baseline dist';
elseif opt == 2
    settings.ylab   = ['Multivariate state transition' newline 'speed (neural dist./ms)'];
    fig_start       = 'speed';
elseif opt == 3
    settings.ylab   = ['Distance between' newline 'trajectories (a.u.)'];
    fig_start       = 'dist bet dur';
end
if is_mult_dur; fig_start = ['Mult-dur ', fig_start]; else; fig_start = [fig_start, ' ci-']; end
settings.fig_name = sprintf('Multivar %s dynamics %s', fig_start, join(settings.ROInames(ROIs)));
end

function signif_times = show_signif(masks, time_vec, ROInames, cat_names)
% masks - cat\dur x time x reg, or time x reg
% prints the time ranges which are significant
if ~exist('cat_names','var'); n_cat = 1; else; n_cat = length(cat_names); end
if ismatrix(masks); masks = permute(masks, [3 1 2]); end
n_reg = size(masks,3);
signif_times = cell(n_cat, n_reg); cat_str = "";
for reg_i = 1:n_reg
for cat_i = 1:n_cat
if n_cat ~= 1; cat_str = cat_names{cat_i} + ": "; end
tmp = bwconncomp(masks(cat_i, :, reg_i)); 
signif_times{cat_i, reg_i} = cell2mat(cellfun(@(x) time_vec(x([1,end])), tmp.PixelIdxList', 'UniformOutput', false));
disp(ROInames{reg_i}+ ": "+cat_str+join(join(string(signif_times{cat_i,reg_i}),"-")'+"ms",", "))
end
end
end

function fig = ci_cat_reg_plot(dat, settings, get_ci_n_val)
% dat: cat x reg x perm
% settings: ROIs (what to plot from dim 2 of data, should be all probably)
%       ylim, ylab, title, cat_cmap, fig_name, ci_pcnt, zero_ln, ROInames_full
[actual_val, low_ci, high_ci] = get_ci_n_val(dat(:,settings.ROIs,:), 3, settings.ci_pcnt);
n_plotted_regs = length(settings.ROIs); xlims = [0.5 n_plotted_regs+0.5];
errbar_inputs = {'o', 'MarkerEdgeColor', 'none', 'Color', 'k', 'MarkerSize', 8, 'CapSize', 2};

fig = figure('Units','Normalized','Position',[0.25 0.25 0.25 0.4], 'Name', settings.fig_name); hold on;
if settings.zero_ln; plot(xlims,[0 0],'k:'); end
for categ = 1:size(dat,1) % plots all categories
    errorbar((1:n_plotted_regs) + (categ-2.5)/8, actual_val(categ,:), actual_val(categ,:)-low_ci(categ,:), high_ci(categ,:)-actual_val(categ,:),...
        errbar_inputs{:}, 'MarkerFaceColor', settings.cat_cmap(categ,:));
end
for reg = 1:(n_plotted_regs-1);plot([reg,reg]+0.5,settings.ylim,'k:');end
xlim(xlims); ylim(settings.ylim); ylabel(settings.ylab); title(settings.title)
set(gca,'XTick',1:n_plotted_regs,'XTickLabels', settings.ROInames(settings.ROIs),'Box','off','FontSize',14);
end