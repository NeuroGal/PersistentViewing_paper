%% Figures for results section 5 (Figure 5 & Supp Figures 8-11)
% generate data with rsa_calc.m

% ! note I will be assuming that rsa_calc ran with z-scored reliability measures 
% (z_by_perm), with the basic statistics being cluster-based (with firstlevel 
% threshold of 1.5 z-score units, as in the manuscript), and additional   
% permutation statistics which will be used for generalization matrices.

% starting with general settings:

clear all;clc;close all
%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));
load([DATA_FOLDER,'figure_settings.mat'])
metric = 'z-score'; % I am assuming this is what was run in the main
get_save_str = @(ctype,dis,subs,esamp) sprintf('%s_%s_subs%s_%s',ctype,dis,join(string(subs),''),esamp);
get_fig_ext = @(dis,subs,esamp) sprintf('%s-%dsub-%s', dis(1:3), length(subs), esamp);

%% Loading the main computation data (5 patient group, correlation, all responsive electrodes)
%   used in: Figure 5d-f (more data will be loaded below for d-e), Figure S8b-c, Figure S9a-b

% Refer to rsa_calc.m to see what these settings mean - do not change them here, we have additional loading below.
subs = [1 2 5 7 8]; e_subsample = 'resp'; dissimilarity = 'correlation';
calc_type      = 'standard';
main_ROIs_load = [1:6,9:10]; % only 1:4 are needed for Figure 5, the rest are for supplementary
main_fig_ext = get_fig_ext(dissimilarity,subs,e_subsample);
general_save_name = get_save_str(calc_type,dissimilarity,subs,e_subsample);

[main_reg_dat,main_settings] = loader(DATA_FOLDER, general_save_name, ROInames, main_ROIs_load);
get_fig_title = @(ROIs,rel_i,cur_tit,fig_ext) sprintf('RSA %s %s-rel %s %s',fig_ext,upper(main_settings.reliability{rel_i}),cur_tit,join(ROInames(ROIs),'-'));
time_vec = main_settings.time_vec;

% plot settings
gen_plot_inputs_line = {'largest_clust','chance',0,'ylab',metric};
diag_viz_dat = []; diag_viz_dat.plot_inputs = gen_plot_inputs_line; diag_viz_dat.font_size = font_size;
diag_viz_dat.ROInames_full = ROInames_full; diag_viz_dat.time_vec = time_vec; diag_viz_dat.tit_cmap = reg_cmap;
gen_plot_inputs_mat  = {'largest_clust','chance_level',0,'metric',metric};
mat_viz_dat = []; mat_viz_dat.plot_inputs = gen_plot_inputs_mat; mat_viz_dat.metric = metric;
mat_viz_dat.ROInames_full = ROInames_full; mat_viz_dat.time_vec = time_vec;

%% Load model data (includes both correlation with model RDMs and reliability with model RDMs partialled)
%   assumed you started with the load above
%   used in: Figure 5d-e, Figure S10

calc_type         = 'model';
mod_ROIs_load     = 1:6; % only 1:4 are needed for Figure 5, the rest are for supplementary
general_save_name = get_save_str(calc_type,dissimilarity,subs,e_subsample);
[mod_reg_dat,mod_settings,model_corr] = loader(DATA_FOLDER, general_save_name, ROInames, mod_ROIs_load);

% some settings
model_rdm_types = mod_settings.model_rdm_types; n_model = length(model_rdm_types);
model_cmap = cat(3,0.2*ones(size(reg_cmap)),reg_cmap_dark,reg_cmap_bright,repmat(cat_cmap(1,:),length(reg_cmap),1)); 
[peak_corr_best,best_mod_idx] = max(max(model_corr.corrs(time_vec>0,:,:),[],1),[],2);
peak_corr_best = squeeze(peak_corr_best); best_mod_idx = squeeze(best_mod_idx); % index of the most correlated model per ROI

% info needed for Figure S10a captions
var_not_exp = 100*(1-model_corr.corrs.^2); tp_thresh = 100;
mean_var_not_exp = squeeze(mean(var_not_exp(time_vec>=tp_thresh,:,:),1))'; % reg x mod
exp_var_best = cellfun(@(row,idx) row(idx), num2cell(mean_var_not_exp,2), num2cell(best_mod_idx));
array2table([string(model_rdm_types(best_mod_idx))',round(peak_corr_best,2),round(exp_var_best)],'VariableNames',{'model','peak','%unexp var'},'RowNames',ROInames(1:max(mod_ROIs_load)))

% create merged diag stats with the original data without partialling anything out added as a sort of first model
%   used in Figure 5d-e and Figure S10b-c
diag_main_n_models = [main_reg_dat.diags(:,:,1:max(mod_ROIs_load));mod_reg_dat.diags(:,:,1:max(mod_ROIs_load))]; % model x reliability x region

%% Figure 5d-e (d: ITEM, e: GEOM)
%   Original (unpartialled) reliability (reg_cmap) and reliability with the most correlated model partialled out (gray)

ROIs = 1:4;
RELs = 1:2;  % what to plot from the reliability metrics (1: ITEM, 2: GEOM)
tmp = arrayfun(@(reg,b_idx) diag_main_n_models([b_idx+1 1],:,reg),(1:max(ROIs))',best_mod_idx(1:max(ROIs)),'UniformOutput',false);
curr_diag_info = cat(3,tmp{:}); 
fig_title = 'dynamics no category info removed vs most correlated model removed';
diag_viz_dat.cmap = cat(3,0.65*ones(size(reg_cmap)),reg_cmap);
get_fig_title_cur = @(rel_i) get_fig_title(ROIs,rel_i,fig_title,main_fig_ext);
figs = rsa_diag_plots(ROIs, RELs, curr_diag_info, diag_viz_dat, get_fig_title_cur);

%% Figure 5f
%   Item reliability generalization in time
%   x\y labels were changed in illustrator

ROIs = 1:2; rel_i = 1;
all_mats = main_reg_dat.relis(:,:,:,rel_i,ROIs);
clims = [-5 10]; % you can also use limer, this is just what I used
tmp = (100:200:899)'; mat_ts = fliplr(num2cell([tmp,tmp+200],2)'); clear tmp
plot_inputs = [gen_plot_inputs_mat,{'no_diag','clims',clims,'time_slices',mat_ts}];
for reg_i = ROIs
    rel_met = all_mats(:,:,:,:,ROIs==reg_i); rel_mask = main_reg_dat.masks_perm(:,:,:,rel_i,reg_i);
    fig = plot_results_mat(rel_met, time_vec, ROInames_full(reg_i), 'mask', rel_mask, plot_inputs{:}, 'line_col', reg_cmap(reg_i,:));
    set(fig,'Name', get_fig_title(reg_i,rel_i,'matrix',main_fig_ext));
end

%% Figure S8b-c
%   Standard (no category info removed) item & Geom reliability.
%   [Figure S8a was a brain w electrodes]
%   b: SM & LT c: OFC & LPFC
%   the next panels will be plotted belows
%   In the manuscript I also plotted marginal clusters (p<0.1), you can
%   redo the stats with rdm_corr_stats.m changing the p_thresh setting.

RELs = 1:2; diag_viz_dat.cmap = reg_cmap; diag_viz_dat.ylims = [-5 10; -3 8]; % this is what I used
for panel = 1:2
    ROIs = 1+4*panel+(0:1);
    get_fig_title_cur = @(rel_i) get_fig_title(ROIs,rel_i,'diagonal dynamics',main_fig_ext);
    figs = rsa_diag_plots(ROIs, RELs, main_reg_dat.diags, diag_viz_dat, get_fig_title_cur);
end

%% Figure S9a
%   Plot mean diagonal & generalization values & % relative to peak

ROIs = 3:6; RELs = 1:2; tp_thresh = 100; 
txt_dat = []; txt_dat.xticks = ROInames(ROIs); txt_dat.ylab = metric; txt_dat.ylab2 = '% from peak'; 
bar_viz_dat = []; bar_viz_dat.ylims2 = [0 100]; bar_viz_dat.fig_w = 0.02; bar_viz_dat.cmap = cat(3, reg_cmap_dark(ROIs,:), reg_cmap_bright(ROIs,:));

for rel_i = RELs
txt_dat.tit = sprintf('%s reliability (> %dms)',main_settings.reliability{rel_i},tp_thresh);
if rel_i == 1; bar_viz_dat.ylims = [0 10]; elseif rel_i == 2; bar_viz_dat.ylims = [0 8]; end % this is what I used
prop_vals = nan(length(ROIs), 2); prop_vals2 = nan(length(ROIs), 2);
for reg_i = ROIs
    cur_peak = max(diag(main_reg_dat.relis(time_vec>0,time_vec>0,:,rel_i,reg_i))); % this is what I used in the paper, similar for other definitions
    thresh_rel = main_reg_dat.relis(time_vec>tp_thresh,time_vec>tp_thresh,:,rel_i,reg_i);   
    prop_vals(ROIs==reg_i,1) = mean(diag(thresh_rel)); % mean_diag
    prop_vals(ROIs==reg_i,2) = mean(thresh_rel(~eye(size(thresh_rel)))); % mean_gen
    prop_vals2(ROIs==reg_i,:) = 100*prop_vals(ROIs==reg_i,:)./cur_peak;
end
fig = bar_plot_prop(prop_vals, [], txt_dat, bar_viz_dat, prop_vals2);
set(fig,'Name', get_fig_title(ROIs,rel_i,sprintf('mean diag gen %dms',tp_thresh),main_fig_ext));
end

%% Figure S9b
%   Geometry reliability of VT & Occ (complementing Figure 5f)

ROIs = [2,1]; rel_i = 2; mat_viz_dat.cmap = reg_cmap; mat_viz_dat.clims = [-5 10; -4 8]; % this is what I used
get_fig_title_cur = @(rel_i) get_fig_title(ROIs,rel_i,'matrices',main_fig_ext);
figs = rsa_mat_plots(ROIs, rel_i, main_reg_dat, mat_viz_dat, get_fig_title_cur);

%% Figures S8d-f, S9c-e: Load additional data followed immediately by plotting
%   Figure S8: Standard item & Geom reliability for ROIs 1:4
%   Figure S9: Item & geometry generalization for ROIs 2:3
%   Relevant settings:
%       Figures S8d,S9c: 5subj,pos,correlation
%       Figures S8e,S9d: 5subj,resp,euclidean
%       Figures S8f,S9e: 8subj,resp,correlation

calc_type = 'standard'; cur_ROIs_load = 1:4; RELs = 1:2;
for setting_i = 1:3
    if setting_i == 1
        cur_subs = [1 2 5 7 8]; cur_esamp = 'pos'; cur_dis = 'correlation';
    elseif setting_i == 2
        cur_subs = [1 2 5 7 8]; cur_esamp = 'resp'; cur_dis = 'euclidean';
    elseif setting_i == 3
        cur_subs = [1 2 5:10]; cur_esamp = 'resp'; cur_dis = 'correlation';
    end
    cur_save_name = get_save_str(calc_type,cur_dis,cur_subs,cur_esamp);
    cur_reg_dat = loader(DATA_FOLDER, cur_save_name, ROInames, cur_ROIs_load);
    cur_fig_ext = get_fig_ext(cur_dis,cur_subs,cur_esamp);

    % Figure S8
    ROIs = 1:4; diag_viz_dat.cmap = reg_cmap; diag_viz_dat.ylims = [-5 10; -3 8]; % this is what I used
    get_fig_title_cur = @(rel_i) get_fig_title(ROIs,rel_i,'diagonal dynamics',cur_fig_ext);
    figs = rsa_diag_plots(ROIs, RELs, cur_reg_dat.diags, diag_viz_dat, get_fig_title_cur);
    
    % Figure S9
    ROIs = [2,1]; mat_viz_dat.cmap = reg_cmap; mat_viz_dat.clims = [-5 10; -4 8]; % this is what I used
    get_fig_title_cur = @(rel_i) get_fig_title(ROIs,rel_i,'matrices',main_fig_ext);
    figs = rsa_mat_plots(ROIs, RELs, cur_reg_dat, mat_viz_dat, get_fig_title_cur);
end

%% Figure S10a
%   data RDM-model RDM correlation dynamics
%   both figures here we used to compose S10a (different y-axis)

for fig_i = 1:2
if fig_i == 1; ROIs = 1:2; else; ROIs = 3:6; end
nrow = 1; ncol = length(ROIs); fig_pos = [0.1 0.25 0.025+0.155*ncol 0.21];
marg_h = 0.135*[2 1]; gap = 0.02; marg_w = 0.01*([5 1]+[2 1]/ncol);
[fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w);
ylims = limer(model_corr.corrs(:,:,ROIs),10); ylims(1) = ylims(1)-0.125*range(ylims);
plot_inputs = {'chance',0,'ylims',ylims,'ylab',['Spearman' newline 'Correlation (\rho)']};
for reg_i = ROIs
    axes(ha(ROIs==reg_i));
    for mod_i = 1:n_model
        corr_r = model_corr.corrs(:,mod_i,reg_i); mask = model_corr.masks(:,mod_i,reg_i); cluster_p = model_corr.clust_ps{mod_i,reg_i};
        nice_line_plot(corr_r, time_vec, '', plot_inputs{:}, 'color', model_cmap(reg_i,:,mod_i)); hold on;
        add_mask_stars(mask, cluster_p, corr_r, time_vec, ylims, model_cmap(reg_i,:,mod_i), [20/range(time_vec) 0.035*mod_i], 'left', 11.5);
    end
    if reg_i ~= ROIs(1); set(gca,'Xlabel',[],'Ylabel',[],'YTick',[]); else; cur_ticks = get(gca,'YTick'); if cur_ticks(end)~=ylims(2); set(gca,'YTick',[get(gca,'YTick'),ylims(2)]); end; end
    title(ROInames_full(reg_i),'Color',reg_cmap(reg_i,:))
end
set(ha,'FontSize',font_size);
set(fig,'Name',sprintf('RSA %s RDM-model correlation dynamics %s',main_fig_ext,join(ROInames(ROIs),'-')));
end

%% Figure S10b-c (b: ITEM, c: GEOM)
%   Original (unpartialled) reliability (gray) with all of the model-partialled results (model_cmap)
%   I changed which plot is on top in illustrator

ROIs = 1:6; RELs = 1:2;  
curr_diag_info = diag_main_n_models;
fig_title = 'dynamics no category info removed vs most correlated model removed';
diag_viz_dat.cmap =cat(3,0.65*ones(size(reg_cmap)),model_cmap); 
get_fig_title_cur = @(rel_i) get_fig_title(ROIs,rel_i,fig_title,main_fig_ext);
figs = rsa_diag_plots(ROIs, RELs, curr_diag_info, diag_viz_dat, get_fig_title_cur);

%% Figure S10d-e (d: ITEM, e: GEOM)
%   Generalization after partialling category models
%   Using the plotting function with models instead of regions
%   In d the x-y labels are accidentally flipped (should be 'time same repetition' in y-label)

ROIs = [2,1]; RELs = 1:2;
cur_viz_dat = mat_viz_dat;
cur_viz_dat.ROInames_full = [model_rdm_types,{'With category information'}]; 
cur_viz_dat.ylims = [-4 10; -4 8]; % optional, this is what I used

for reg_i = ROIs
cur_viz_dat.cmap = [squeeze(model_cmap(reg_i,:,:))';0.65*ones(size(reg_cmap))]; % shift so models is the first dimension
cur_mod_dat = []; % structure so models in the place of regions (and add the standard computation in the end)
cur_mod_dat.relis = permute(cat(3,mod_reg_dat.relis(:,:,:,:,reg_i),main_reg_dat.relis(:,:,:,:,reg_i)),[1 2 5 4 3]);
cur_mod_dat.masks_perm = permute(cat(3,mod_reg_dat.masks_perm(:,:,:,:,reg_i),main_reg_dat.masks_perm(:,:,:,:,reg_i)),[1 2 5 4 3]);
get_fig_title_cur = @(rel_i) get_fig_title(reg_i,rel_i,'model matrices',main_fig_ext);
figs = rsa_mat_plots(1:(n_model+1), RELs, cur_mod_dat, cur_viz_dat, get_fig_title_cur);
end

%% Load single-category results (Figure S11)

calc_type         = 'categories';
cat_ROIs_load     = 1:4;
general_save_name = get_save_str(calc_type,dissimilarity,subs,e_subsample);
[cat_reg_dat,cat_settings] = loader(DATA_FOLDER, general_save_name, ROInames, cat_ROIs_load);
cat_names = cat_settings.cat_names; n_cat = length(cat_names);
cat_rel_i = 1;

% settings for Figure S11a-c
ex_cat_ROIs = [1 2 4]; chosen_cat = [4 3 2 1]; peak_times = [450 255 0 195]; % 0 is just a filler, this can be easily computed so it's not really needed to hardcode it
cat_rng = 1:max(ex_cat_ROIs);
ex_cat_titles = cellfun(@(str1,str2) [str1,' - ',upper(str2(1)),str2(2:end),'s'], cellstr(ROInames_full(cat_rng)), cat_names(chosen_cat), 'UniformOutput',false);
ex_cat_mergerA = @(field) permute(arrayfun(@(reg,categ) cat_reg_dat.(field){categ,cat_rel_i,reg},cat_rng, chosen_cat(cat_rng),'UniformOutput',false),[1 3 2]);
ex_cat_mergerB = @(field) cell2mat(permute(arrayfun(@(reg,categ) cat_reg_dat.(field)(:,:,categ,cat_rel_i,reg), cat_rng, chosen_cat(cat_rng),'UniformOutput',false),[1 3 4 5 2]));

%% Figure S11a
%   State-space visualization for one example category in each region.
%   Uses code from state_space_plot repository. See  plot_mv_patterns for more options.
%   Uses distinguishable_colors https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors (added to helper functions folder)

add_images = true; % this is saved on our lab server, so you need to plot it without
plot_settings  = []; plot_settings.single = true; plot_settings.std = false; plot_settings.use_sem = false; plot_settings.mean = false;
viz_settings   = []; viz_settings.new_fig = false; viz_settings.show_ticks = false;
image_settings = []; image_settings.folder = 'S:\Lab-Shared\Experiments\Duration_Gamma\ECOG\Original Exp\exp\stimuli\';
cat_tmp = {'face','watch','',''}'; all_image_filenames = string(cat_tmp(cat_settings.stim_info.categ)) + pad(string(cat_settings.stim_info.stim_id),3,'left','0') + ".bmp";

for reg_i = ex_cat_ROIs
cat_i = chosen_cat(reg_i); peak_t = peak_times(reg_i); stim_idx = cat_settings.stim_info.categ == cat_i; n_stim = sum(stim_idx)/2;
image_settings.filenames = all_image_filenames(stim_idx);
plot_patt_inputs  = {'plot_settings',plot_settings,'viz_settings',viz_settings};
if add_images; plot_patt_inputs = [plot_patt_inputs,{'images',image_settings}]; end
fig_title = sprintf('MDS tsne %s %s %dms', ROInames(reg_i), cat_settings.cat_names{cat_i}, peak_t);
fig = figure('Name',fig_title,'Units','normalized','Position',[0.05 0.05 0.35 0.5]);
ha = axes('Units','normalized','Position',[0.05 0.05 0.9 0.8]);
rng(5) % for reproducability
dim_red_patterns = do_dim_red(cat_settings.segs(stim_idx,cat_settings.ROI_idx==reg_i,:),'tsne',find(time_vec == peak_t),'metric', dissimilarity,'n_dim', 2);
plot_mv_patterns(dim_red_patterns, peak_t, 'color_code',grp2idx(categorical(cat_settings.stim_info.stim_id(stim_idx))),'colormap',distinguishable_colors(n_stim), plot_patt_inputs{:});
title(ha,ex_cat_titles(reg_i),'Visible','on','Units','Normalized','Position',[0.5 1.12 0]);
set(ha,'Visible','off')
end

%% Figure S11b
%   One example category dynamics (Occ-animals,VT-objects,PFC-faces)
%   Recolor title in illustrator & add vertical horizonal line

curr_diag_info = ex_cat_mergerA('diags');
cur_viz_dat = diag_viz_dat; cur_viz_dat.ROInames_full = ex_cat_titles;
cur_viz_dat.cmap = reg_cmap; cur_viz_dat.ylims = [-3.5 5.5]; % optional
get_fig_title_cur = @(rel_i) get_fig_title(ex_cat_ROIs,rel_i,'example single category dynamics',main_fig_ext);
fig = rsa_diag_plots(ex_cat_ROIs, cat_rel_i, curr_diag_info, cur_viz_dat, get_fig_title_cur);
fig.Children(end).YTick = -2.5:2.5:cur_viz_dat.ylims(2);

%% Figure S11c
%   One example category matrices (Occ-animals,VT-objects,PFC-faces)
%   color-axis was a bit different in the MS, you can tweak that

curr_reg_dat = []; curr_reg_dat.relis = ex_cat_mergerB('relis');
curr_reg_dat.masks = ex_cat_mergerB('masks'); curr_reg_dat.clust_ps = ex_cat_mergerA('clust_ps');
cur_viz_dat = mat_viz_dat; cur_viz_dat.ROInames_full = ex_cat_titles;
cur_viz_dat.cmap = reg_cmap;
get_fig_title_cur = @(rel_i) get_fig_title(ex_cat_ROIs,rel_i,'example single category matrices',main_fig_ext); use_clust = true;
fig = rsa_mat_plots(ex_cat_ROIs, cat_rel_i, curr_reg_dat, cur_viz_dat, get_fig_title_cur, use_clust);

%% Figure S11d
%   region x category dynamics

ROIs = 1:4;
all_diags = cellfun_wrap(@(x) x.vals,permute(cat_reg_dat.diags(:,cat_rel_i,ROIs),[3 1 2]),true);
ylims = limer(all_diags,1);
ncol = n_cat; nrow = length(ROIs);
gap = [0.015+0.04/length(ROIs) 0.02]; marg_h =[0.04 0.02]+[0.48 0.24]/length(ROIs); marg_w = [0.22 0.02]; fig_size = [0.05 0.1 0.525 0.065+0.075*length(ROIs)];
[fig,ha] = fig_ax_wrap(fig_size, nrow, ncol, gap, marg_h, marg_w);
for reg_i = ROIs
    for categ = 1:n_cat
        p_idx = (reg_i-1)*ncol + categ; cur_diag = cat_reg_dat.diags{categ,cat_rel_i,reg_i};
        axes(ha(p_idx));
        rel_met = cur_diag.vals; mask = cur_diag.mask; cluster_p = cur_diag.cluster_p;
        plot_inputs = [gen_plot_inputs_line,{'ylims',ylims,'color',reg_cmap(reg_i,:),'mask',mask,'cluster_pval',cluster_p}];
        nice_line_plot(rel_met, time_vec, '', plot_inputs{:}); hold on
        if reg_i ~= ROIs(end); set(gca,'XTick',[]); end
        if p_idx ~= ncol*(nrow-1)+1; set(gca,'XLabel',[],'YLabel',[],'YTick',[]); end
        if reg_i == 1; title(gca,cat_names{categ},'Color',cat_cmap(categ,:)); end
        if categ == 2
            ylb = ylabel(ROInames_full(reg_i),'Color',reg_cmap(reg_i,:));
            set(ylb,'Rotation',0,'FontWeight','bold','HorizontalAlignment','center')
            ylb.Position(1) = -2000;
        end
    end
end
set(fig,'Name',get_fig_title(ex_cat_ROIs,cat_rel_i,'regions x categories dynamics',main_fig_ext));

%% functions

function [reg_dat,settings,model_corr] = loader(DATA_FOLDER, save_name, ROInames, ROIs_load)
% information for all regions
settings = load([DATA_FOLDER, sprintf('rsa_settings_%s.mat', save_name)]);
n_time = length(settings.time_vec); n_reg = max(ROIs_load); n_rel = length(settings.reliability);
is_model = false; is_categ = false;
switch save_name(1:5)
    case 'stand'
        dim = 1;
    case 'model'
        dim = length(settings.model_rdm_types); % n_model
        is_model = true;
    case 'categ'
        dim = length(settings.cat_names); % n_cat
        is_categ = true;
end
% load region specific data - saving in the regions axis will be according to the constant order
dims = [n_time, n_time, dim, n_rel, n_reg]; relis = nan(dims); masks = nan(dims); masks_perm = nan(dims);
dims = [dim, n_rel, n_reg]; clust_ps = cell(dims); diags = cell(dims);
model_corr = []; 
if is_model % currently not saving pt by pt perms since it's not needed
    dims = [n_time, dim, n_reg];
    model_corr.corrs = nan(dims); model_corr.masks = nan(dims); model_corr.clust_ps = cell(dim,n_reg);
end
for reg_i = ROIs_load
    load([DATA_FOLDER,sprintf('rsa_%s_%s.mat',ROInames(reg_i),save_name)])
    relis(:,:,:,:,reg_i) = cell2mat(permute(corr_r(:,1:n_rel), [3 4 1 2])); % only to n_rel since the last entry is correlation with the model and the dimensions are different
    masks(:,:,:,:,reg_i) = cell2mat(permute(mask(:,1:n_rel), [3 4 1 2])); 
    if ~is_categ; masks_perm(:,:,:,:,reg_i) = cell2mat(permute(mask_perm(:,1:n_rel), [3 4 1 2])); end
    clust_ps(:,:,reg_i) = cluster_p(:,1:n_rel);
    diags(:,:,reg_i) = diag_stats; % here there is no entry for the model correlations so it's fine
    if is_model
        model_corr.corrs(:,:,reg_i) = cell2mat(corr_r(:,n_rel+1)');
        model_corr.masks(:,:,reg_i) = cell2mat(mask(:,n_rel+1)');
        model_corr.clust_ps(:,reg_i) = cluster_p(:,n_rel+1);
    end
    fprintf('Loaded %s (#%d/%d)\n',ROInames(reg_i),find(ROIs_load==reg_i),length(ROIs_load))
end
reg_dat = []; fields = {'relis','masks','masks_perm','clust_ps','diags'};
for f = 1:length(fields); eval(sprintf('reg_dat.%s = %s;', fields{f}, fields{f})); end
end