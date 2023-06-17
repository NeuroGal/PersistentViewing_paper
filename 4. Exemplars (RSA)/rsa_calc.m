%% RSA\single exemplar analyses: results section 5 (Figure 5 & Supp Figures 8-11)
% take the output from here to rsa_fig.m
clear all;clc;close all

%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));

%% settings

% the specific_segs.m code has a randomization stage to decide which
% repetition will be first\second, to take the exact segments I used in the
% paper instead of segments you generated with that use this:
paper_segs    = true;
% ! for the single categories you have to use this & set smooth to 100
% (specific segs has data smoothed with 50 ms moving average not 100 ms).
smooth        = 50; % in ms, moving average smoothing (50\100)

% subject group: -> defines which segments we load
%   Use [1,2,5,7,8] for almost everything - 5subj
%   [1:2,5:10] for Figure S8f and Figure S9e - 8subj
subs          = [1 2 5 7 8];

% which electrodes? 'resp'(responsive)\'pos' (responsive & increasing from baseline) [can also use 'all', but it wasn't in the paper]
%   Use 'resp' for almost everything
%   'pos' for Figure S8d and Figure S9c
e_subsample   = 'resp';

% what dissimilarity metric? 'correlation','euclidean','mahalanobis','cosine' - whatever goes into matlab's pdist (only the first two were used in the paper)
%   Use 'correlation' for almost everything
%   'euclidean' for Figure S8e and Figure S9d
dissimilarity = 'correlation';

% in practice these are the only combinations you need for these settings:
%   Figures S8d,S9c: 5subj,pos,correlation | Figures S8e,S9d: 5subj,resp,euclidean | Figures S8f,S9e: 8subj,resp,correlation
%   everything else: 5subj,resp,correlation

% what calculation to run? ('standard'\'model'\'categories')
%   settings: 5subj,resp,correlation
%       'standard'   - correlating the RDM of repetition 1 with repetition 2 in the same & other time points
%                      ROIs=[1:6,9:10] runs with both reliability measures (Figures 5d-f,S9a-b require 1:4; Figures S9a,S8b-c require 5:6,9:10)
%       'model'      - correlating with category model RDMs (see Figure 5c), and repeating the standard calculation after partialling out the model data
%                      ROIs=1:6 runs with both reliability measures (Figures 5d-e (only 1:4 needed), S10 (1:6 needed))
%                      the reliability metric only matters for the partialling out, the correlation is done with a full averaged RDM in either case so this is done only once
%       'categories' - run the standard calculation in each category separately
%                      ROIs=1:4 runs with 'item' only (Figure S11)
%   for Figures S8d-f, S9c-e run 'standard',ROIs=1:4,both reliability measures with the subject,electrode&dissimilarity settings explained above
calc_type     = 'standard';

% which ROIs - order: Occ,VT,Par,PFC,SM,LT,Occ-retin,Occ-not-retin,OFC,LPFC
%   probably 1:4\1:6\1:10 (or [1:6,9:10]), see details under what calculation to run
ROIs          = [1:6,9:10];

%% load data (created with specific_segs.m):

load([DATA_FOLDER,'figure_settings.mat'],'ROInames')
% constant settings
n_rep = 2; d_str = 'long'; downsample_hz = 200; cat_nums = 1:4; av_in_cat = false;
segs_desc = sprintf('%drep_dur%s_subs%s_%dHz_cat%s_av%s',n_rep,d_str,join(string(subs),''),...
    downsample_hz,join(string(cat_nums),''),string(av_in_cat));
load_name = sprintf('segs_%s.mat',segs_desc);
if paper_segs; load_name = sprintf('paper%dms_%s',smooth,load_name); end
load([DATA_FOLDER, load_name]);
switch e_subsample
    case 'resp'
        e_filt = elec_info.IsResponsive;
    case 'all'
        e_filt = true(size(elec_info.IsResponsive));
    case 'pos'
        e_filt = elec_info.IsResponsive & elec_info.RespSign==1;
end
segs = segs(:,e_filt,:); elec_info = elec_info(e_filt,:);

general_save_name = sprintf('%s_%s_subs%s_%s',calc_type,dissimilarity,join(string(subs),''),e_subsample);

%% constant calculation settings

corr_type = 'Spearman';  %  'Pearson'\'Kendall'\'Spearman' (the default in the function) [doesn't matter much, the results were similar, I used Spearman in the paper]
z_by_perm = true; % z-score the rdm correlations using the permutation distribution (this is how I defined the reliability metrics in the paper) - doesn't apply to the model correlated with the rdm

model_rdm_types = {'Single-category', 'Low-level', 'Semantic', 'Face-vs-rest'}; % see Figure 5c
n_models = length(model_rdm_types);

stat_info = [];
% stat_type: 'analytic'\'perm'\'cluster'\'max', how to do correction for multiple comparisons
%   cluster uses cluster-based permutations, max with max-statistic control for MC the other two are corrected with Benjamini-Hochberg FDR
%   note I added the max option only after completing the manuscript so it wasn't tested thoroughly
%   don't use 'analytic' with z_by_perm
stat_info.stat_type      = 'cluster'; % I use cluster here and compute FDR when needed in the figures code for matrices etc.
stat_info.n_sides        = 1; % !!! right sided (larger than)\dual sided test
stat_info.p_thresh       = 0.05; % threshold for significance (after fdr\for the clusters)
stat_info.n_perm         = 1000;
stat_info.get_diag_stats = true;
stat_info.chance_level   = 0; % true for correlation & z-score
% cluster specific settings (see inside the function)
stat_info.n_clusters        = Inf;
stat_info.firstlevel_type   = 'stat'; % actually z-score if z-by-perm
stat_info.firstlevel_thresh = 1.5;    % in z-score units for this case, note this is changed below for the categories calc_type to account for the added smoothing

% specific for correlating with model RDM
stat_info_model_corr = stat_info;
stat_info_model_corr.firstlevel_type   = 'p';
stat_info_model_corr.firstlevel_thresh = 0.05;
stat_info_model_corr.get_diag_stats    = false;

% we repeat the statistics with point by point permutations here which was
% used for some of the figures and then we don't need to save the permutations themsevles.
stat_info_perm = stat_info;
stat_info_perm.stat_type = 'perm';

inputs_corr_calc = {'corr_type', corr_type};
if z_by_perm; inputs_corr_calc = [inputs_corr_calc, {'z_by_perm'}]; end

% reliability measure (see the manuscript for details, Figure 5a-b): 'item','geom'
general_save = {'time_vec','reliability'};
switch calc_type
    case 'standard'
        reliability   = {'item','geom'};
    case 'model' % currently not saving the models themselves, nor stat_info_model_corr
        reliability   = {'item','geom'};
        general_save  = [general_save,{'model_rdm_types'}];
    case 'categories'
        reliability   = {'item'}; % you can add 'geom', I just didn't do it in the paper (and currently there is no loop over reliability in that section below)
        general_save  = [general_save,{'cat_names','segs','stim_info','ROI_idx'}]; % segs,stim_info & ROI_idx added for S11a (assuming max(ROIs) was 6, if higher this needs to be adjusted) 
        stat_info.firstlevel_thresh = stat_info.firstlevel_thresh/2;
        stat_info.p_thresh          = 0.06; % ! I show a marginal cluster in Figure S11, you can change this
end
n_rel = length(reliability); n_cat = length(cat_names);

%% calculation (takes time, less than decoding)

for reg_i = ROIs
if reg_i<=6; ROI_idx = grp2idx(elec_info.ROI); else; ROI_idx = grp2idx(elec_info.ROISplit); end % use ROIsplit if any of regions 7-10 requested
segs_roi = segs(:,ROI_idx==reg_i,:); n_elec = size(segs_roi,2);
if n_elec <= 2 ; continue; end % nothing to run on 

rdm = calculate_rdm(segs_roi, dissimilarity, false); 
n_time = size(rdm,3);n_stim = length(categories_vec)/2;
tmp = mat2cell(rdm,[n_stim n_stim],[n_stim n_stim],n_time); 
mean_rdm = mean(cat(4, tmp{:}),4);

switch calc_type
    case 'standard'
        dims = [1,n_rel];
        corr_r = cell(dims); mask = cell(dims); cluster_p = cell(dims); diag_stats = cell(dims); mask_perm = cell(dims); % currently saved
        corr_p = cell(dims); corr_perm_r = cell(dims); % currently not saved
        for rel_i = 1:n_rel
        rng(10) % for reproducability (not really crucial, that's just how I ran it for the paper, without it numbers are very similar but not identical)
        [corr_r{rel_i}, mask{rel_i}, corr_p{rel_i}, cluster_p{rel_i}, corr_perm_r{rel_i}, diag_stats{rel_i}] = ...
            rdm_rel_calc(rdm, reliability{rel_i}, stat_info, inputs_corr_calc{:});
        [mask_perm{rel_i}, ~, ~, ~, diag_stats_tmp] = rdm_corr_stats(stat_info_perm, corr_r{rel_i}, corr_perm_r{rel_i}); % I'm not saving the pvalues, FDR threshold etc, add if needed for some reason
        diag_stats{rel_i}.mask_perm = diag_stats_tmp.mask;
        end
    case 'model'
        dims = [n_models,n_rel+1];
        corr_r = cell(dims); mask = cell(dims); cluster_p = cell(dims); mask_perm = cell(dims); % currently saved
        diag_stats = cell(n_models,n_rel); % also saved - no diag for the corr (it's a single line anyway)
        corr_p = cell(dims); corr_perm_r = cell(dims); % currently not saved
        for model_i = 1:n_models 
            model_rdm = get_model_rdm(categories_vec(1:n_stim), cat_names, model_rdm_types{model_i});
            % compute rdm correlations partialling the model
            for rel_i = 1:n_rel
            rng(10) % for reproducability
            [corr_r{model_i,rel_i}, mask{model_i,rel_i}, corr_p{model_i,rel_i}, cluster_p{model_i,rel_i}, corr_perm_r{model_i,rel_i}, diag_stats{model_i,rel_i}] = ...
                rdm_rel_calc(rdm, reliability{rel_i}, stat_info, inputs_corr_calc{:}, 'model_rdm', model_rdm);
            [mask_perm{model_i,rel_i}, ~, ~, ~, diag_stats_tmp] = rdm_corr_stats(stat_info_perm, corr_r{model_i,rel_i}, corr_perm_r{model_i,rel_i}); % I'm not saving the pvalues, FDR threshold etc, add if needed for some reason
            diag_stats{model_i,rel_i}.mask_perm = diag_stats_tmp.mask;
            end
            % compute correlation with model
            rng(10) % for reproducability - originally I put the rng(10) once before a loop with all correlations with models so the numbers will be a bit different
            [corr_r{model_i,3}, mask{model_i,3}, corr_p{model_i,3}, cluster_p{model_i,3}, corr_perm_r{model_i,3}] = corr_w_cat_mod(mean_rdm, model_rdm, stat_info_model_corr, corr_type);
            % not using point-by-point permutations anywhere here so not saving that.
        end
    case 'categories'
        rng(10) % for reproducability (originally I had it before all categories so leaving that here for now)
        dims = [n_cat,1]; % assuming only item reliability is used
        corr_r = cell(dims); mask = cell(dims); cluster_p = cell(dims); mask_perm = cell(dims); diag_stats = cell(n_models,n_rel); % currently saved
        corr_p = cell(dims); corr_perm_r = cell(dims); % currently not saved
        for categ = [1,3,4,2] % using this instead of 1:n_cat because it was my original order, but it really doesn't matter (only for me to get the exactly same values with the clean saving structure)
            this_cat_rdm = rdm(categories_vec==categ,categories_vec==categ,:);
            if numel(this_cat_rdm)>1 % make sure it's not a single dist\ empty
                [corr_r{categ}, mask{categ}, corr_p{categ}, cluster_p{categ}, corr_perm_r{categ}, diag_stats{categ}] = ...
                    rdm_rel_calc(this_cat_rdm, reliability{1}, stat_info, inputs_corr_calc{:});   
            end
        end
        % no need to save point-by-point permutation data it's not used in the manuscript
end
save([DATA_FOLDER,sprintf('rsa_%s_%s.mat',ROInames(reg_i),general_save_name)],'corr_r','mask','cluster_p','mask_perm','diag_stats') 
% this is pretty minimalistic - just what we need for rsa_fig.m, you can add corr_p, corr_r_perm, rdm, mean_rdm etc. for completeness.
end

% saving settings relevant for all ROIs (to be used later with the figures) - also very minimalistic, 
% add stat_info, z_by_perm, categories_vec etc if you want it to be complete.
save([DATA_FOLDER, sprintf('rsa_settings_%s.mat', general_save_name)],general_save{:});