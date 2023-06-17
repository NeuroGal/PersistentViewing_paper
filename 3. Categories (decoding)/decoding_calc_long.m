%% Decoding analyses for all long (>=900ms) stimuli:
% Result section 3 (Figure 3 & Supp Figures 4-6)
% Take the output from here to decoding_fig_long.m

clear all;clc;close all

%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));

%% settings

% specific controls (overrides setting below)
subsample_control = false; % true only for reproducing Figure S4c (electrode subsampling for VT and Occ)
bipolar_control   = false; % true only for reproducing Figure S6d (bipolar-unipolar reference comparison for OFC)
% for Figure S6e (ofc with\without saccades decoding comparison) see ofc_saccade_control.m

% Define what segments we load (from specific_segs.m - run that before):
% Subjects:
%   most figures (Figure 3, Figure S4, S5c-e, S6c): subs = 1:10
%   single patients (S5a-b, S6a-b): change subs to a single patient each time (i.e. subs = 1, then subs = 2 etc.)
%   special controls: 
%       Figure S6d: subs = [8, 10]
%       Figure S6e: subs = 8 (see ofc_saccade_control.m)
subs        = 1:10;

% Categories: (face,watch,object,animal)
%   almost everything: cat_nums = 1:4 - this will be used in decoding_wrapper.m to do 1vs1 decoding for all pairs
%                                       (the function can also do multiclass with the right input)
%   Figure 5a-b: cat_nums = 1:2 - just run the face-watch comparison and save yourself a lot of time
cat_nums    = 1:4;

% Define what electrodes\ROIs we use to do the decoding:
%   (all non-noisy electrodes & all ROIs are loaded in the initial segments)

% which electrodes? 'all' (all electrodes)\'resp' (responsive)\'pos' (responsive & increasing from baseline)
%   Use 'resp' for almost everything
%   'all' & 'pos' are needed for Figure S4a, 'all' for Figure S6d-e
e_subsample = 'resp';

% ROIs: (Occ,VT,Par,PFC,SM,LT,Occ-retin,Occ-not-retin,OFC,LPFC)
% which are used where:
%   Figure 3: 1:4
%   Figure S4: a) 1:6, b) 1:4, c) 1:2
%   Figure S5: a,d) 2 b,e) 1 c) 7:8
%   Figure S6: a) 3 b) 4 c) 9:10 d-e) specific electrodes in region 4 (PFC)
% in practice:
%   with all patients - run 1:10 for e_sample = 'resp' (used in main & supp), 1:6 for e_sample = 'all'\'pos' (just used for Figure S4)
%   single patients: ROIs = 1:2 for Figure S5a-b, and ROIs = 3:4 for Figure S6a-b
%   PFC controls: specific electrodes in ROI = 4, see below and ofc_saccade_control.m
ROIs        = 1:10;

if subsample_control; subs = 1:10; cat_nums = 1:4; e_subsample = 'resp'; ROIs = 1:2; n_subsample = 36; n_iter = 1000; end % 36 = number of PFC electrodes
if bipolar_control; subs = [8, 10]; cat_nums = 1:4; e_subsample = 'all'; ROIs = [4 4]; run_bipolar = true; end % ROIs = 4 is written twice, once it will run with unipolar referencing and once with bipolar referencing

%% load data (created with specific_segs.m):

load([DATA_FOLDER,'figure_settings.mat'],'ROInames')
% constant settings
d_str = 'long'; downsample_hz = 200; n_rep = 1; av_in_cat = false;
segs_desc = sprintf('%drep_dur%s_subs%s_%dHz_cat%s_av%s',n_rep,d_str,join(string(subs),''),...
    downsample_hz,join(string(cat_nums),''),string(av_in_cat));
load([DATA_FOLDER, sprintf('segs_%s.mat',segs_desc)]);

% some more followups on the settings
switch e_subsample
    case 'resp'
        e_filt = elec_info.IsResponsive;
    case 'all'
        e_filt = true(size(elec_info.IsResponsive));
    case 'pos'
        e_filt = elec_info.IsResponsive & elec_info.RespSign==1;
end
segs = segs(:,e_filt,:); elec_info = elec_info(e_filt,:);
n_cat = length(cat_nums);

general_save_name = sprintf('dur%s_subs%s_%s_cat%s',d_str,join(string(subs),''),e_subsample,join(string(cat_nums),''));

%% decoding settings
%   Decoding uses MVPA_Light, and recieves the settings structure used
%   there (so see that for documentation)
%       Treder, MS Frontiers in Neuroscience (2020) https://doi.org/10.3389/fnins.2020.00289
%       Code here: https://github.com/treder/MVPA-Light
%   * stats are different! uses custom written code for single-subject (in
%   our case mega-patient) stats.

% decoding settings: see MVPA Light
cfg_decoding = [];
cfg_decoding.preprocess = 'undersample';
cfg_decoding.metric     = 'auc';

% statistics settings: see decoding_stats.m for more details
cfg_stats = [];
% stat_type: 'analytic'\'perm'\'cluster'\'max', how to do correction for multiple comparisons
%   cluster uses cluster-based permutations, max with max-statistic control for MC
%   the other two are corrected with Benjamini-Hochberg FDR
%   note I added the max option only after completing the manuscript so it wasn't tested thoroughly
cfg_stats.stat_type    = 'cluster'; % I use cluster here and compute FDR when needed in the figures code for matrices etc.
cfg_stats.n_sides      = 1; % !!! right sided (larger than)\dual sided test
cfg_stats.p_thresh     = 0.05; % threshold for significance (after fdr\for the clusters)
cfg_stats.n_perm       = 1000;
cfg_stats.acc          = strcmp(cfg_decoding.metric,'acc'); % this may be overwritten inside decoding_wrapper.m
cfg_stats.chance_level = 0.5; % specific for our case since we're running 2 balanced classes one vs  the other

% cluster specific settings (see inside the function)
cfg_stats.n_clusters        = Inf;
cfg_stats.firstlevel_type   = 'stat';
cfg_stats.firstlevel_thresh = 0.1; % note that the thresh is relative to chance level!

% we repeat the statistics with point by point permutations here which was
% used for some of the figures and then we don't need to save the permutations themsevles
% (we do need their diagonals to compare between regions for Figure S4B, but it's not as heavy)
cfg_stats_perm           = cfg_stats;
cfg_stats_perm.stat_type = 'perm';

%% calculation (can take some time)
% The loop over category pairs happens within the function decoding_wrapper.m
% For each comparison we compute the temporal generalization matrix (King &
%   Dehaene, 2014), and the 'standard' decoding is the diagonal of this.
%   stats are computed for both so we won't 'over-correct' for multiple
%   comparisons if we don't use the full temporal generalization matrix.

for reg_i = ROIs
if reg_i<=6; ROI_idx = grp2idx(elec_info.ROI); else; ROI_idx = grp2idx(elec_info.ROISplit); end % use ROIsplit if any of regions 7-10 requested
segs_roi = segs(:,ROI_idx==reg_i,:); n_elec = size(segs_roi,2); save_name = ROInames(reg_i);
if n_elec < 1; continue; end % nothing to run on (e.g. 1:3 PFC there are no responsive electrodes)

if bipolar_control % all of this is just for running the control for figure s6d
    cur_elec_idx = [8:18,23:34]'; % indices of the ofc strips excluding one bad electrode 
    segs_roi = segs_roi(:,cur_elec_idx,:); % ofc strips data
    orig_elec_idx = [49:58,60,25:36]'; % electrode indices numbered as they were within each patient (s8 until 60 then s10)
    adj_idx = cellfun_wrap(@(x) intersect(x,orig_elec_idx),num2cell([orig_elec_idx+1,orig_elec_idx-1,orig_elec_idx+6,orig_elec_idx-6],2)); % adjacent electrodes on the strips
    new_idx = cellfun_wrap(@(e_idx) find(ismember(orig_elec_idx,e_idx)),adj_idx); % move to the numbers of the current segs_roi
    if run_bipolar
        segs_roi = segs_roi - cellfun_wrap(@(x) mean(segs_roi(:,x,:),2), new_idx', true);
        save_name = 'bipolar'; run_bipolar = false; % so the next time will be unipolar
    else
        save_name = 'unipolar';
    end
    cur_elec_resp = elec_info.IsResponsive(ROI_idx==reg_i,:); cur_elec_resp = cur_elec_resp(cur_elec_idx);
    segs_roi = segs_roi(:, cur_elec_resp, :); % only done here since we wanted to rereference also with the non-responsive neighbors
end

if n_cat>2
% >>> just to match the paper results 1:1 <<<
% for some reason when I ran it the order of the categories in the segments was faces,objects,animals,watches (not faces,watches,objects,animals)
% I am changing this here to get exactly the values in the figures (with the random seed), but it's very similar without it
% I will flip it back after the computation of the cross-validation & permutations
categories_vec_tmp = grp2idx(categorical(categories_vec,[1,3,4,2]));
[categories_vec_tmp,tmp_idx]=sort(categories_vec_tmp);
segs_roi_tmp = segs_roi(tmp_idx,:,:); reorder_idx = [3 1 2 5 6 4];
else
categories_vec_tmp = categories_vec; segs_roi_tmp = segs_roi;
end

if ~subsample_control
rng(10) % for reproducability
[perf_all, stat_all, cat_order, perms_all] = decoding_wrapper(segs_roi_tmp, categories_vec_tmp, cfg_decoding, cfg_stats); % the actual magic happens in the function

if n_cat>2
% >>> reorder results to the original order of faces,objects,animals,watches - 
% the comparison order above is: F-O,F-A,F-W,O-A,O-W,A-W. We flip it to: F-W,F-O,F-A,W-O,W-A,O-A
perf_all = perf_all(:,:,reorder_idx); stat_all = stat_all(:,:,reorder_idx,:); perms_all = perms_all(:,reorder_idx);
end

% fixes for how I presented the data in the paper
perf_all    = squeeze(100*perf_all); % turn it to % success (not decimal out of 1) -> time x time x n_comparisons
perms_all   = cellfun_wrap(@(x) 100*x, perms_all); % 1 x n_comparisons -> saving here only the diagonals to use later, can easily save the full permutations it's just heavy
% extracting stats to be used in the decoding_fig_long.m:
clust_masks = squeeze(stat_all(1,:,:,:)); % (!) I am assuming the basic statistics above are run with cluster permutations, cell: comparisons x matrix\diag
cluster_p   = squeeze(stat_all(2,:,:,:)); % same
if n_cat == 2; clust_masks = clust_masks'; cluster_p = cluster_p'; end
% adding point-by-point permutations corrected for MC with FDR (for the figures, avoid saving the full permutations like this):
perm_masks  = cell(size(clust_masks));
for comp_i = 1:size(perf_all,3)
perm_masks{comp_i,1} = decoding_stats(cfg_stats_perm, perf_all(:,:,comp_i), perms_all{comp_i});
perm_masks{comp_i,2} = decoding_stats(cfg_stats_perm, diag(perf_all(:,:,comp_i)), get_diags(perms_all{comp_i}));
end
perms_diag = cellfun_wrap(@(x) get_diags(x), perms_all, true); % time x comparison x permutation
save([DATA_FOLDER,sprintf('decoding_%s_%s.mat',save_name,general_save_name)],'perf_all', 'clust_masks', 'cluster_p', 'perm_masks', 'perms_diag')
% saved just what we need for the figures, you can also add segs_roi, full permutations, FDR threshold etc for completeness
else % yes subsample_control
save_name = sprintf('%s%delec',save_name,n_subsample);
rng(10) % for reproducability
e_idx_all = false(n_elec ,n_iter); % electrodes used in each repetition
for rep = 1:n_iter; e_idx_all(randsample(n_elec,n_subsample), rep) = true; end
rng(10) % for reproducability
perf_all = nan(length(time_vec), length(time_vec), nchoosek(length(cat_names),2), n_iter);
for rep = 1:n_iter
    perf_all(:,:,:,rep) = decoding_wrapper(segs_roi_tmp(:,e_idx_all(:,rep),:), categories_vec_tmp, cfg_decoding); % no stats\permutations needed here so no cfg_stats
end
perf_all = 100*permute(get_diags(perf_all(:,:,reorder_idx,:)),[1 4 3 2]); % time x repetition x comparison, % success
save([DATA_FOLDER,sprintf('decoding_%s_%s.mat',save_name,general_save_name)],'perf_all', 'e_idx_all') % we really just need perf_all
end
end

if ~subsample_control
% saving settings relevant for all ROIs (to be used later with the figures)
cfg_stats.firstlevel_thresh = 10; % adjust to this units change
cfg_stats.chance_level = 50; % same
if n_cat > 2
cat_order = nchoosek(1:4,2); % this should have also been the output above, it's changed only if some categories are skipped because they are too small.
% Since we changed the order back after the change to categories_vec_tmp we can use the original cat_names and this will create the correct component names
else
cat_ord = [1 2];
end
comp_names      = join(cellfun_wrap(@(x) [upper(x(1)),x(2:end)],cat_names(cat_order)),'-'); % cat_order is the same for all regions so that's why we only save it here
save([DATA_FOLDER, sprintf('decoding_settings_%s.mat', general_save_name)],'comp_names','time_vec','cfg_stats');
% you can add cfg_decoding & categories_vec, they are just not used in the figure generating code.
end