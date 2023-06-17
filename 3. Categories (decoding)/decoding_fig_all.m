%% Figures for results section 4 (Figure 4 & Supp Figure 7)

clear all;clc;close all
%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));
load([DATA_FOLDER,'figure_settings.mat'])

subs 		= 4:10; % 1:3 or 4:10
% 4:10 used in Figure 4 (regions occ,vt,par - see below) and Supp Fig 7c-f (same regions)
% 1:3  used in Figure 4a-b and Supp Fig 7c-f for Par and Supp Fig 7b-e for Occ,VT

% mostly stable settings
e_subsample = 'resp'; % all panels except Supp Figure 7e which uses 'pos'
cat_nums    = 1:2;    % face-watch, all panels except Supp Figure 7d which uses [1,3] (face-object, animals did not have enough exemplars)

% constant settings
d_str = 'all';
if all(ismember(subs,4:10)); ROIs_load = [1:2,4]; else; ROIs_load = 1:3; end % order is: occ,vt,par,pfc 
general_save_name = sprintf('dur%s_subs%s_%s_cat%s',d_str,join(string(subs),''),e_subsample,join(string(cat_nums),''));

%% Loading data

% load general settings
load([DATA_FOLDER, sprintf('decoding_settings_%s.mat', general_save_name)]);
chance_res = 50; chance_diff = 0;

% load ROI specific data
data_all   = cell(1,max(ROIs_load)); % saving in the regions axis will be according to the constant order
for reg_i = ROIs_load
    load([DATA_FOLDER,sprintf('decoding_%s_%s.mat',ROInames(reg_i),general_save_name)])
    data_all{reg_i} = data_out;
    fprintf('Loaded %s (#%d/%d)\n',ROInames(reg_i),find(ROIs_load==reg_i),length(ROIs_load))
end
clear data_out 

%% Diagonals plot - Fig 4A-B, Supp Fig 7B-E

isperm = false; % false for everything except Figure 7c

all_diags = permute(cellfun_wrap(@(x) get_diags(x.res_mats), data_all(ROIs_load), true),[1 3 2]); % time x dur x reg
all_diffs = diff(all_diags,[],2);
ylims_diags = limer(all_diags, 1/5); ylims_diffs = limer(all_diffs, 1/5);

fig = figure('Units','Normalized','Position',[0.1 0.25 0.54 0.32]);
nrow = 1; ncol = length(ROIs_load); marg_h = [0.19 0.105]; marg_w = [0.055 0.015]; gap = 0.0125; gap_bet_ha = 0.035; height = (1-sum(marg_h))/2;
ha_diag = tight_subplot(nrow, ncol, gap, [marg_h(1)+gap_bet_ha+height marg_h(2)], marg_w);
ha_diff = tight_subplot(nrow, ncol, gap, [marg_h(1) marg_h(2)+gap_bet_ha+height], marg_w);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]

diag_cols = cat(3, adj_cmap(reg_cmap, 1), adj_cmap(reg_cmap, 0.8), adj_cmap(reg_cmap,0.6));
diff_cols = cat(3, adj_cmap(reg_cmap, 0.9), adj_cmap(reg_cmap, 0.7));
ln_dist = 0.04;
fig_inputs_diag = {'ylab','AUC (%)','chance',chance_res,'ylims',ylims_diags}; 
fig_inputs_diff = {'ylab','AUC diff (%)','chance',chance_diff,'ylims',ylims_diffs};

for reg_i = ROIs_load
    rel_dat = data_all{reg_i};
    for plt = 1:2 % diag, diff
        if plt == 1 % diag
            vals = all_diags(:, :, ROIs_load==reg_i); pvals = rel_dat.pvals_res(:,2);
            if ~isperm; mask = rel_dat.res_diag_masks; else; mask = rel_dat.res_diag_masks_perm; end
            ax = ha_diag; ylims = ylims_diags; cur_cols = diag_cols; dur_i_idx = [3,2,1];
            fig_inputs = fig_inputs_diag; ln_offset_func = @(dur_i) ln_dist*(4-dur_i); 
        else % diff
            vals = all_diffs(:, :, ROIs_load==reg_i); pvals = rel_dat.pvals_clust2;
            if ~isperm; mask = rel_dat.diff_diag_masks_clust2; else; mask = rel_dat.diff_diag_masks_perm; end
            ax = ha_diff; ylims = ylims_diffs; cur_cols = diff_cols; dur_i_idx = 1:2;
            fig_inputs = fig_inputs_diff; ln_offset_func = @(dur_i) ln_dist*(3-dur_i);
        end
        axes(ax(ROIs_load==reg_i))    
        for dur_i = 1:3; plot([600*(dur_i-0.5) 600*(dur_i-0.5)], ylims, 'Color', diag_cols(reg_i,:,dur_i)); hold on; end
        for dur_i = dur_i_idx
            if ~isperm; cluster_pvals = pvals{dur_i}; else; cluster_pvals = []; end % it's not really cluster_pvals in this case
            col = cur_cols(reg_i,:,dur_i); offset = [1, ln_offset_func(dur_i)];
			nice_line_plot(vals(:,dur_i), time_vec, '', fig_inputs{:}, 'color',col)
            add_mask_stars(mask(:,:,dur_i), cluster_pvals, vals(:, dur_i), time_vec, ylims, col, offset, 'right', 11.5)
        end
		if plt == 1; title(ROInames_full(reg_i),'FontSize',14,'Color',reg_cmap(reg_i,:)); end
    end
end
set([ha_diag(2:end) ha_diff(2:end)],'YTick',[],'YLabel',[],'XLabel',[],'XTick',[]);
set(ha_diff(1),'XTick',[0,300:600:2100],'YTick', [ylims_diffs(1), 0, ylims_diffs(2)]);
set(ha_diag(1),'XLabel',[],'XTick', [], 'YTick', [ylims_diags(1), 50, ylims_diags(2)])
set(fig,'Name',sprintf('Decoding multiple durations dynamics and differences %s subs-%s isperm-%s', comp_name, join(string(subs),''), string(isperm)));

%% Matrix plot - Fig 4C-D, Supp Fig 7F

if all(ismember(subs,4:10)); ROIs_plot = [1,2,4]; else; ROIs_plot = 3; end
for reg_i = ROIs_plot
dat = data_all{reg_i}; res_mats = dat.res_mats; diff_mats = diff(res_mats,[],3);
cluster_pvals = {dat.pvals_res(:,1), dat.pvals_diff(:,1)};
suptitle = []; suptitle.txt = ROInames_full(reg_i); suptitle.col = reg_cmap(reg_i,:);
fig = mult_dur_mat_plot(res_mats, diff_mats, time_vec, dat.res_masks, dat.diff_masks, cluster_pvals, suptitle);
set(fig,'Name',sprintf('Decoding multiple durations temporal generalization matrices %s',ROInames(reg_i)));
end

%% functions

function cmap = adj_cmap(cmap, bright_const)
cmap_hsv = rgb2hsv(cmap);
cmap_hsv(:,3) = bright_const;
cmap = hsv2rgb(cmap_hsv);
end