%% Figures for results section 3 (Figure 3 & Supp Figures 4-6)
% also includes final statistics for some of the figures

clear all;clc;close all
%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));

load([DATA_FOLDER,'figure_settings.mat']); d_str = 'long'; 
get_save_str = @(subs,esamp,cats) sprintf('durlong_subs%s_%s_cat%s',join(string(subs),''),esamp,join(string(cats),''));
fig_namer = @(fig, title_txt, ROIs) set(fig,'Name',sprintf('Decoding %s %s',title_txt,join(ROInames(ROIs),'-')));
chance_level = 50; % specific for our case since we're running 2 balanced classes one vs the other
metric = 'AUC (%)'; % assuming this is the metric in the main
fig_inputs = [];
fig_inputs.line = {'largest_clust','chance',chance_level,'ylab',metric};
fig_inputs.mat  = {'decoding','chance_level',chance_level,'metric',metric};

%% Loading the main computation data (all patients, responsive electrodes) - load for settings for other figures too
%   used in: Figure 3, Figure S4a-b, Figure S5c-e, Figure S6c
%   more data for Figure S4a,c will be loaded below
%   generate with decoding_calc_long.m

ROIs_load = 1:10; % only 1:4 are needed for Figure 3, the rest are for supplementary
subs = 1:10; e_subsample = 'resp'; cat_nums = 1:4; 
main_save_name = get_save_str(subs,e_subsample,cat_nums);

[main_reg_dat,main_settings] = loader(DATA_FOLDER, main_save_name, ROInames, ROIs_load);
time_vec = main_settings.time_vec; comp_names = main_settings.comp_names; n_comp = length(comp_names); cfg_stats = main_settings.cfg_stats;

%% Main Figure
%% Figure 3B: peaks across comparisons
% uses the permutation values for computing stats here

ROIs = 1:4; 
peak_vals = cellfun(@(x) max(x,[],1), main_reg_dat.results(:,ROIs,2)); % ncomp x 1
perm_peak = cellfun_wrap(@(x) max(x,[],1), main_reg_dat.diag_perms(:,ROIs), true); % ncomp x 1 x nperm
pvals_peak  = (sum(peak_vals<perm_peak,3)+1)/(size(perm_peak,3)+1);
thresh_peak = prctile(perm_peak, 95, 3);
stat_dat = []; stat_dat.pvals = pvals_peak; stat_dat.thresh = thresh_peak; stat_dat.add_stars = true; stat_dat.add_lines = true;
viz_dat = []; viz_dat.cmap = permute(reg_cmap(ROIs,:),[3 2 1]); viz_dat.ylims = [50 105]; viz_dat.fig_w = 0.018;
txt_dat = []; txt_dat.tit = 'Peak decoding'; txt_dat.ylab = 'AUC (%)'; txt_dat.xticks = comp_names; txt_dat.add_legend = true; txt_dat.leg = ROInames_full(ROIs);
fig = bar_plot_prop(peak_vals, stat_dat, txt_dat, viz_dat);
fig_namer(fig, 'peaks across comparisons', ROIs);

% calculations for main text
array2table(round([mean(peak_vals,1);std(peak_vals,[],1)/sqrt(n_comp)],1),'VariableNames',ROInames(ROIs),'RowNames',{'mean','sem'})
array2table(pvals_peak,'VariableNames',ROInames(ROIs),'RowNames',comp_names)

%% Figure 3C: face-watch decoding dynamics

comp_i = 1; ROIs = 1:4; 
nrow = 1; ncol = ceil(length(ROIs)/nrow); fig_pos = [0.1 0.25 0.025+0.155*ncol 0.08+0.13*nrow];
marg_h = [0.05+0.22/nrow 0.015+0.12/nrow]; marg_w = [0.05+0.08/ncol 0.01+0.03/ncol]; gap = [0.105+0.08/nrow 0.01];
[fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w);
fig_inputs_cur = [{'ylims',[30 100]},fig_inputs.line];
for reg_i = ROIs
    axes(ha(ROIs==reg_i));
    vals = main_reg_dat.results{comp_i,reg_i,2};
    mask_clust = main_reg_dat.masks_clust{comp_i,reg_i,2}; mask_perm = main_reg_dat.masks_perm{comp_i,reg_i,2};
    pval = main_reg_dat.clust_ps{comp_i,reg_i,2}; % assumes the calculation code already ran cluster based statistics
    nice_line_plot(vals, time_vec, '', fig_inputs_cur{:},'mask',mask_clust,'color',reg_cmap(reg_i,:),'cluster_pval',pval); hold on
    set(findobj(gca,'Type','Text'),'Color','r')    
    nice_line_plot(vals, time_vec, ROInames_full(reg_i), fig_inputs_cur{:}, 'mask_info', struct('color','k','offset',0.035),...
        'mask',mask_perm,'color',reg_cmap(reg_i,:)); % add point by point w FDR lines
    if reg_i ~= ROIs(1)
        set(gca,'XTick',[],'YTick',[],'Xlabel',[],'Ylabel',[]);
    else
        set(gca,'YTick',[50 100]);
    end
end
set(ha,'FontSize',font_size); fig_namer(fig, sprintf('dynamics %s',comp_names{comp_i}), ROIs);

% numbers for the text:
round(cellfun(@(x) mean(x(time_vec>100)), main_reg_dat.results(comp_i,1:2,2)),1)

%% Figure 3D: face-watch temporal generalization matrices (TGMs)

comp_i = 1; ROIs = 1:4; 
tmp = (100:200:899)'; ts = fliplr(num2cell([tmp,tmp+200],2)');
fig_inputs_cur = [fig_inputs.mat,{'no_diag','time_slices',ts,'clims',[0 100]}];
for reg_i = ROIs
    dat = main_reg_dat.results{comp_i,reg_i,1}; mask = main_reg_dat.masks_perm{comp_i,reg_i,1};
    fig = plot_results_mat(dat, time_vec, ROInames_full(reg_i), 'mask', mask, 'line_col', reg_cmap(reg_i,:), fig_inputs_cur{:});
    fig_namer(fig, sprintf('TGM %s',comp_names{comp_i}), reg_i);
end

%% Supplementary
%% Loading for Figure S4(A,C)

% Figure S4A: load results with all electrodes & positive only electrodes
ROIs_load = 1:6;
for e_subsample = ["pos","all"]
cur_save_name = get_save_str(subs,e_subsample,cat_nums);
cur_reg_dat = loader(DATA_FOLDER, cur_save_name, ROInames, ROIs_load);
eval(sprintf('%s_reg_dat = cur_reg_dat;',e_subsample))
end

% Figure S4C: load VT and Occ results after subsampling the number of electrodes
ROIs_subsample = 1:2; n_elec_red = 36; 
results_red = cell(max(ROIs_subsample),1);
for reg_i = ROIs_subsample
    load([DATA_FOLDER,sprintf('decoding_%s%delec_%s.mat',ROInames(reg_i),n_elec_red,main_save_name)])
    results_red{reg_i} = perf_all; % time x repetition x comparison
    fprintf('Loaded %s\n',ROInames(reg_i))
end
clear perf_all

%% Figure S4A: all regions x all comparisons 3 types of decoding: responsive, pos-only, all elec overlayed

ROIs = 1:6; ylims = [0 100]; 
[fig,ha] = fig4_figure(0.55,length(ROIs),n_comp,[0.09 0.05],0.013);
cmap_all = nan([size(reg_cmap),3]); bright_idx = [0.8 1 0.6];
for plt = 1:3; cmap_hsv = rgb2hsv(reg_cmap); cmap_hsv(:,3) = bright_idx(plt); cmap_all(:,:,plt) = hsv2rgb(cmap_hsv); end

for comp_i = 1:n_comp
    for reg_i = ROIs
        plt_idx = (find(ROIs==reg_i)-1)*n_comp+comp_i;
        axes(ha(plt_idx));
        for idx = [1,3,2] % 1- resp, 2- only pos, 3- all
            if idx == 1
                cur_reg_dat = main_reg_dat;
            elseif idx == 2
                cur_reg_dat = pos_reg_dat;
            else
                cur_reg_dat = all_reg_dat;
            end
            vals = cur_reg_dat.results{comp_i,reg_i,2}; mask = cur_reg_dat.masks_clust{comp_i,reg_i,2}; pval = cur_reg_dat.clust_ps{comp_i,reg_i,2};
            col = cmap_all(reg_i,:,idx);
            nice_line_plot(vals, time_vec, '', fig_inputs.line{:},'color',col,'ylims',ylims); hold on
            add_mask_stars(mask, pval, vals, time_vec, ylims, col, [0 0.06*idx],'left')
        end 
        fig4_final_touch(reg_i,comp_i,ROIs,comp_names,ROInames_full,reg_cmap);
    end
end
set(ha,'FontSize',13); fig_namer(fig, 'all regions x all comparisons - responsive, pos-only, all elec', ROIs);

%% Figure S4B: four main regions difference plots

ROIs = [2 1 4 3]; % no good reason, but that's how I plotted it
diff_pairs = nchoosek(ROIs,2);
all_diags = cell2mat(permute(main_reg_dat.results(:,:,2),[3 1 2])); % time x comp x reg  
all_diffs = all_diags(:,:,diff_pairs(:,1))-all_diags(:,:,diff_pairs(:,2));  % time x comp x reg (diff)
ylims = limer(all_diffs, 1/25);
cfg_stats_diff = cfg_stats; cfg_stats_diff.chance_level = 0; cfg_stats_diff.n_sides = 2;
[fig,ha] = fig4_figure(0.55,size(diff_pairs,1),n_comp,[0.09 0.05],0.013);
fig_inputs_cur = [fig_inputs.line,{'ylims',ylims}];
fig_inputs_cur{find(strcmp(fig_inputs.line,'chance'))+1}=0;
for pair_i = 1:size(diff_pairs,1)
    vals = all_diffs(:, :, pair_i);
    perms = diff(cell2mat(permute(main_reg_dat.diag_perms(:, diff_pairs(pair_i,:)),[3 1 4 2])),[],4); % time x comp x perm
    col = mean(reg_cmap(diff_pairs(pair_i,:),:));
    for comp_i = 1:n_comp
        axes(ha(n_comp*(pair_i-1)+comp_i))
        [mask, cluster_pvals] = decoding_stats(cfg_stats_diff, vals(:, comp_i), perms(:, comp_i, :));
        nice_line_plot(vals(:, comp_i), time_vec, '',fig_inputs_cur{:},'color',col)
        add_mask_stars(mask, cluster_pvals, vals(:, comp_i), time_vec, ylims, col, [0 0.06],'left')
        if pair_i == 1; title(comp_names(comp_i),'Color','k'); end
        if comp_i ~= 1 || pair_i ~= size(diff_pairs,1); set(gca,'XTick',[],'YTick',[],'XLabel',[],'YLabel',[]); else; set(gca,'YTick',[-50 0 50]); end
        if comp_i == 2                    
            ylb = ylabel(join(ROInames_full(diff_pairs(pair_i,:)),[newline,'vs',newline])); % recolored later in illustrator
            set(ylb,'Rotation',0,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
            ylb.Position(1) = -1800;
        end
    end
end
set(ha,'FontSize',13); fig_namer(fig, 'difference between regions (all comps)', ROIs);

%% Figure S4C: VT and Occ subsampling electrodes to the number in PFC

[fig,ha] = fig4_figure(0.225,length(ROIs_subsample),n_comp,2.15*[0.1 0.05],0.01);
for reg_i = ROIs_subsample
    dat = results_red{reg_i}; col = reg_cmap(reg_i,:);        
    fig_inputs_cur = [{'ylims',[20 100],'color',col},fig_inputs.line];
    for comp_i = 1:n_comp
        plt_idx = (find(ROIs_subsample==reg_i)-1)*n_comp+comp_i;
        axes(ha(plt_idx));
        nice_line_plot(mean(dat(:,:,comp_i),2), time_vec, '', fig_inputs_cur{:}); hold on
        varplot(time_vec, dat(:,:,comp_i), 'std', 'Color',col)
        fig4_final_touch(reg_i,comp_i,ROIs_subsample,comp_names,ROInames_full,reg_cmap)
    end
end
set(ha,'FontSize',13); fig_namer(fig, sprintf('%d electrodes subsampling',n_elec_red), ROIs_subsample);

%% Load single patient data for Figure S5A-B and S6A-B

e_subsample = 'resp'; ROIs_load = 1:4;
for reg_i = ROIs_load
    if reg_i <= 2; cat_nums=1:2; elseif reg_i <=4; cat_nums=1:4; end % assuming that's what was run
    dims = [nchoosek(length(cat_nums),2), length(subs), 2]; sub_dat = [];
    sub_dat.results = cell(dims); sub_dat.masks_clust = cell(dims);
    sub_dat.masks_perm = cell(dims); sub_dat.clust_ps = cell(dims);
    for subj = subs
    cur_save_name = get_save_str(subj,e_subsample,cat_nums);
    one_sub_dat = loader(DATA_FOLDER, cur_save_name, ROInames, reg_i);
    if ~isempty(one_sub_dat.diag_perms{1,reg_i}) % indicating something was loaded
        sub_dat.results(:,subj,:) = one_sub_dat.results(:,reg_i,:);
        sub_dat.masks_clust(:,subj,:) = one_sub_dat.masks_clust(:,reg_i,:);
        sub_dat.masks_perm(:,subj,:) = one_sub_dat.masks_perm(:,reg_i,:);
        sub_dat.clust_ps(:,subj,:) = one_sub_dat.clust_ps(:,reg_i,:);
    end
    end
    sub_dat.rel_subs = find(cellfun(@(x) ~isempty(x), sub_dat.results(:,:,1)));
    eval(sprintf('%s_sub_dat = sub_dat;',ROInames(reg_i)));
end
clear sub_dat one_sub_dat

%% Supp 5A-B: Single patients VT, Occ TGM & diagonal plot

ROIs = 1:2;
comp_info = []; comp_info.ord = 1; fig_info = []; fig_info.nrow = 3; fig_info.opt_set = []; 
for reg_i = ROIs
eval(sprintf('cur_sub_dat = %s_sub_dat;',ROInames(reg_i)));
if reg_i == 1; fig_info.del_idx = 6; else; fig_info.del_idx = [8 12]; end
fig_info.ncol = ceil(length(cur_sub_dat.rel_subs)/fig_info.nrow);
reg_info = []; reg_info.ord = cur_sub_dat.rel_subs; reg_info.names = "S"+string(cur_sub_dat.rel_subs); reg_info.reg_cmap = reg_cmap(reg_i,:); 
fig = dec_mat_diag_plot(cur_sub_dat, time_vec, reg_info, comp_info, fig_info, fig_inputs);
fig_namer(fig,'single patients TGMs and diagonals',reg_i)
end

%% Supp 6A-B: Single patients PFC, Par dynamics (diagonal) only

ROIs = 3:4;
comp_info = []; comp_info.ord = [1:3,5]; comp_info.names = comp_names;
reg_info = []; subj_mode = true;
for reg_i = ROIs
eval(sprintf('cur_sub_dat = %s_sub_dat;',ROInames(reg_i)));
reg_info.reg_cmap = reg_cmap(reg_i,:); 
if reg_i == 4; reg_info.ord = [8,10]; else; reg_info.ord = [2:5,8:10]; end
fig = dec_diag_plot(cur_sub_dat, time_vec, reg_info, comp_info, subj_mode, fig_inputs.line);
fig_namer(fig,'single patients diagonals',reg_i)
end

%% Supp Fig 5C-E: occipital retinotopic & non-retinotopic split, Occ and VT all comparisons - diagonals & TGMs
% In the paper the diagonals were removed from panels d-e since they are already presented in Supp Fig 4A.

ROIs = [7 8 2 1];
fig_info = []; fig_info.nrow = length(ROIs); fig_info.ncol = n_comp;
opt_set = []; opt_set.right = 1.5; opt_set.w_const = 0.01; opt_set.ratios = [0.3 0.02 0.15]; fig_info.opt_set = opt_set;
reg_info = []; reg_info.ord = ROIs; reg_info.names = ROInames_full; reg_info.reg_cmap = reg_cmap;
comp_info = []; comp_info.ord = 1:6; comp_info.names = comp_names;
fig = dec_mat_diag_plot(main_reg_dat, time_vec, reg_info, comp_info, fig_info, fig_inputs);
fig_namer(fig, 'all comp matrix and diagonal', ROIs);

% numbers for the legend
ROIs = [7 8];
for meas = 1:3 % peak, mean diag>100, mean gen>100
    if meas == 1
        all_comp_dat = cellfun(@(x) max(x), main_reg_dat.results(:,ROIs,2)); txt = 'Peak';
    elseif meas == 2
        all_comp_dat = cellfun(@(x) mean(x(time_vec>100)), main_reg_dat.results(:,ROIs,2)); txt = 'Mean diag after 100ms';
    elseif meas == 3
        tmp = cellfun_wrap(@(x) x(time_vec>100,time_vec>100), main_reg_dat.results(:,ROIs,1));
        all_comp_dat = cellfun(@(x) mean(x(~eye(size(x)))), tmp); txt = 'Mean generalization after 100ms';
    end
    dat = string(round([mean(all_comp_dat);std(all_comp_dat)/sqrt(n_comp)],1))';
    fprintf([txt,' >>>\n']); disp(ROInames(ROIs)'+": "+dat(:,1)+ " +- "+dat(:,2))
end

%% Supp Fig 6C: plot Par and PFC all comparisons only diagonals (no TGM)

ROIs = 9:10;
reg_info = []; reg_info.ord = ROIs; reg_info.names = ROInames_full; reg_info.reg_cmap = reg_cmap;
comp_info = []; comp_info.ord = [1:3,5]; comp_info.names = comp_names; subj_mode = false;
fig = dec_diag_plot(main_reg_dat, time_vec, reg_info, comp_info, subj_mode, fig_inputs.line);
fig_namer(fig, 'all comp diagonal', ROIs);

% numbers for the legend
all_comp_dat = cellfun(@(x) max(x), main_reg_dat.results(:,ROIs,2));
dat = string(round([mean(all_comp_dat);std(all_comp_dat)/sqrt(n_comp)],1))';
fprintf('Peak >>>\n'); disp(ROInames(ROIs)'+": "+dat(:,1)+ " +- "+dat(:,2))
    
%% Supp Fig 6D: bipolar vs average reference comparison (OFC strips)
% load
% settings
cat_nums = 1:4; subs = [8 10]; e_subsample = 'all'; % not really 'all', but that's what we used for bipolar referencing before taking out only those centered on responsive contacts
fig_name = ['decoding_%s_' get_save_str(subs,e_subsample,cat_nums) '.mat'];
leg_txt = {'Average','Bipolar'};
% data
unipolar = load([DATA_FOLDER, sprintf(fig_name, 'unipolar')]);
bipolar = load([DATA_FOLDER, sprintf(fig_name, 'bipolar')]);

% plot
comp_nums = 1:5;
all_diags = cellfun_wrap(@(x) get_diags(x.perf_all(:,:,comp_nums)),{unipolar,bipolar}, true); % time x reference x comp
ylims = limer(all_diags, 1/5); ylims(2) = 100; ylims(1) = ylims(1)-10;
all_r_cor = cellfun(@(x) corr(x(:,1),x(:,2)), squeeze(num2cell(cellfun_wrap(@(x) get_diags(x.perf_all),{unipolar,bipolar}, true),[1 2])));
[fig, ha] = one_row_fig(length(comp_nums));
for comp_i = comp_nums
    axes(ha(comp_nums==comp_i));
    for line_idx = 1:2 % unipolar, bipolar
        if line_idx == 1; dat = unipolar; else; dat = bipolar; end
        vals = diag(dat.perf_all(:,:,comp_i)); mask = dat.clust_masks{comp_i,2}; pval = dat.cluster_p{comp_i,2};
        nice_line_plot(vals, time_vec, '', fig_inputs.line{:}, 'ylims',ylims,'color',ofc_cols(line_idx,:));
        add_mask_stars(mask, pval, vals, time_vec, ylims, ofc_cols(line_idx,:), [20/range(time_vec) 0.06*line_idx-0.02], 'left', 10)
    end
    title(comp_names{comp_i},'Color','k')
    text(time_vec(end),ylims(2),sprintf('r = %0.2g',all_r_cor(comp_i)),'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top');
    if comp_i ~= comp_nums(1); set(gca,'XTick',[],'YTick',[],'XLabel',[],'YLabel',[]); end
    tmp = findobj(gca,'Type','Line');
    line_handles = tmp(arrayfun(@(x) ~any(isnan(x.YData)), tmp) & arrayfun(@(x) x.LineStyle=='-', tmp));
    if comp_i == comp_nums(end); leg = legend(flipud(line_handles), leg_txt,'Box','off','NumColumns',2); leg.Position(2) = 0.06; end
end
fig_namer(fig, 'comparison bipolar-average reference', ROInames=="OFC");
    
% numbers for legend
[mean(all_r_cor),std(all_r_cor)/n_comp]
array2table(all_r_cor','VariableNames',comp_names)

%% Supp Fig 6E: comparison of decoding in all trials vs trials without any saccades (per time-window)

% load
load([DATA_FOLDER,'decoding_ofc_saccade_control.mat'])

% I: decoding bars
comp_nums = 1:5;
[fig, ha] = one_row_fig(length(comp_nums));
lims = [10 90];
for comp_i = comp_nums
    axes(ha(comp_nums==comp_i)); b_handles = [];
    for b = 1:2
        b_handle(b) = bar((0.4:3.4)+0.4*b,squeeze(perf_all(comp_i,:,b)), 0.35,'LineStyle','none','FaceColor','flat');hold on
        b_handle(b).CData = ofc_cols(b,:);
        add_txt = @(idx, sig, col) text(b_handle(b).XEndPoints(idx),(lims(2)-5)*ones(1,sum(idx)),...
                repmat('*',1,sig),'FontWeight','bold','HorizontalAlignment','center','FontSize',12,'Color',col);
        for sig = 1:3
            signif = pvals_all(comp_i,:,b)<=signif_levels(sig)& pvals_all(comp_i,:,b)>signif_levels(sig+1);
            add_txt(signif & pvals_all(comp_i,:,b) <= 0.05/4, sig, 'red') % bonferroni correction
            add_txt(signif & pvals_all(comp_i,:,b) > 0.05/4, sig, 'black')
        end
    end
    plot([0 5],[chance_level chance_level],'k:'); hold on
    ylim(lims); xlim([0.4 4.6]); title(comp_names{comp_i})
    set(gca,'Box','off','FontSize',font_size,'XTick',[],'YTick',[])
    if comp_i == comp_nums(1); set(gca,'XTick',1:4,'XTickLabel',win_names,'YTick',[lims(1), chance_level, lims(2)]);xlabel('Time-window (ms)');ylabel(metric); end
    if comp_i == comp_nums(end); leg = legend(b_handle,{'All trials', 'No saccades'},'Box','off','NumColumns',2); leg.Position(2) = 0.06; end
end
fig_namer(fig, 'bars comparison all-trials vs no saccades', ROInames=="OFC");

% II: number of trials of each type
ncol = 4; nrow = 1; fig_pos = [0.1 0.2 0.55 0.24]; gap = 0.02; marg_h =  [0.3 0.1]; marg_w = [0.1 0.2];
[fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w);

for w_idx = 1:4
    axes(ha(w_idx))
    b_handle = bar([-1,1:4],squeeze(trial_counts(:,w_idx,[1 1:4]))','stacked','CData',[ones(1,3)*132/255;cat_cmap],... % 1 twice for the color of the legend
        'LineWidth',2,'FaceColor','flat','EdgeColor','flat');
    set(b_handle(1),'FaceColor',[1 1 1]); set(b_handle(2),'FaceAlpha',0.5);

    title(win_names(w_idx)); xlim([0.2 4.8])
    set(gca,'Box','Off','XTick',1:4,'XTickLabels',cat_names,'FontSize',font_size)
    if w_idx == 1; ylabel('# unique images'); else; set(gca,'YTick',[],'XTick',[]); end
    if w_idx == 4; leg = legend({'No saccades','Saccades in some trials','Saccades in all trials'},'Box','off'); leg.Position(1) = 0.78; end
end
fig_namer(fig, 'saccade control number of unique images', ROInames=="OFC");
% legend text
pcnt_excluded = squeeze(100*trial_counts(3,:,:)./sum(trial_counts,1)); pcnt_excluded = pcnt_excluded(:);
round([min(pcnt_excluded) max(pcnt_excluded) mean(pcnt_excluded) std(pcnt_excluded)/sqrt(numel(pcnt_excluded))])

% III: diff all-no sac
%   Jitter is not identical to the paper
fig = figure('Units','Normalized','Position',[0.2 0.2 0.1 0.325]);
diff_vals = -diff(perf_all,[],3); diff_vals = diff_vals(:); 
diff_mean = mean(diff_vals); diff_sem = std(diff_vals)/sqrt(length(diff_vals));
scatter(0.4*rand(size(diff_vals))+0.8, diff_vals, 55, repmat(reg_cmap(9,:),length(diff_vals),1), 'filled', 'MarkerFaceAlpha',0.6); hold on
errorbar(1, diff_mean, diff_sem,'Color','k','LineWidth',1.25);
plot([0.5 1.5],[diff_mean diff_mean],'Color','k','LineWidth',2.25)
set(gca,'XTick',[],'Box','off','FontSize',13.5)
ylim(limer(diff_vals, 1/20));ylabel('All minus no-saccades'); xlim([0 2]);
fig_namer(fig, 'difference (dots) all-trials vs no saccades', ROInames=="OFC");

% IV: corr plot (all vs no sac)
var1 = perf_all(:,:,1); var2 = perf_all(:,:,2); rel_pvals = pvals_all(:,:,2);
fig = figure('Position',[100 100 370 370]);
desc = repmat(comp_names,1,4); lims = [30 90];
plot(lims, lims,'k:'); hold on; plot(lims, [chance_level chance_level],'k:'); plot([chance_level chance_level], lims,'k:'); 
add_scat = @(idx, col) scatter(var1(idx), var2(idx), 50, col, 'filled');
add_scat(rel_pvals>0.05, reg_cmap(9,:)); add_scat(rel_pvals<=0.05/4, 'r')
add_scat(rel_pvals<=0.05 & ~(pvals_all(:,:,2)<=0.05/4), 'k')
idx = rel_pvals<=0.05/4; text(var1(idx)+2,var2(idx), desc(idx))
xlabel('All trials (AUC %)'); ylabel('No saccades (AUC %)')
xlim(lims); ylim(lims)
set(gca,'Box','off','FontSize',font_size,'XTick',lims(1):20:lims(2), 'YTick',lims(1):20:lims(2))
fig_namer(fig, 'correlation all-trials vs no saccades', ROInames=="OFC");
% legend text
[r, p]=corr(var1(:),var2(:))

%% functions

function [fig,ha] = one_row_fig(ncol)
nrow = 1; fig_pos = [0.05 0.05 0.1*(1+ncol) 0.08*(1+nrow)];
gap = [0.008 0.005]; marg_h =  [0.15 0.1]; marg_w = [0.2 0.05]; marg_w = marg_w/(sum(marg_w)*(1+ncol)); marg_h = marg_h/(sum(marg_h)*(1+nrow));
[fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w);
end

function [reg_dat,settings] = loader(DATA_FOLDER, save_name, ROInames, ROIs_load)
% information for all regions
settings = load([DATA_FOLDER, sprintf('decoding_settings_%s.mat', save_name)]);
n_comp = length(settings.comp_names); n_reg = max(ROIs_load);

% load region specific data - saving in the regions axis will be according to the constant order
dims = [n_comp, n_reg, 2]; % last dim: full TGM, diag
results = cell(dims); masks_clust = cell(dims); masks_perm = cell(dims); clust_ps = cell(dims); diag_perms = cell(n_comp,n_reg); % no matrix perms saved
for reg_i = ROIs_load
    FILE_PATH = [DATA_FOLDER,sprintf('decoding_%s_%s.mat',ROInames(reg_i),save_name)]; rel_file = dir(FILE_PATH);
    if isempty(rel_file)
    fprintf('%s data not found\n', ROInames(reg_i));
    else
    load(FILE_PATH)
    results(:,reg_i,1) = num2cell(perf_all,[1 2]);
    results(:,reg_i,2) = cellfun_wrap(@diag,results(:,reg_i,1));
    masks_clust(:,reg_i,:) = clust_masks;
    masks_perm(:,reg_i,:)  = perm_masks;
    clust_ps(:,reg_i,:)    = cluster_p;
    diag_perms(:,reg_i)    = num2cell(perms_diag,[1 3]);
    fprintf('Loaded %s (#%d/%d)\n',ROInames(reg_i),find(ROIs_load==reg_i),length(ROIs_load))
    end
end
reg_dat = []; fields = {'results','masks_clust','masks_perm','clust_ps','diag_perms'};
for f = 1:length(fields); eval(sprintf('reg_dat.%s = %s;', fields{f}, fields{f})); end
end

function fig4_final_touch(reg_i,comp_i,ROIs,comp_names,ROInames_full,reg_cmap)
if reg_i == ROIs(1); title(comp_names(comp_i),'Color','k'); end
if comp_i ~= 1 || reg_i ~= ROIs(end)
    set(gca,'XTick',[],'YTick',[],'XLabel',[],'YLabel',[]);   
else
    set(gca,'YTick',[0 50 100]); 
end
if comp_i == 2
    ylb = ylabel(ROInames_full(reg_i),'Color',reg_cmap(reg_i,:));
    set(ylb,'Rotation',0,'FontWeight','bold','HorizontalAlignment','center')
    ylb.Position(1) = -1800;
end
end

function [fig,ha] = fig4_figure(fig_height,nrow,ncol,marg_h,gap1)
gap = [gap1 0.009]; marg_w = [0.15 0.01]; fig_size = [0.05 0.1 0.675 fig_height];
[fig,ha] = fig_ax_wrap(fig_size, nrow, ncol, gap, marg_h, marg_w);
end