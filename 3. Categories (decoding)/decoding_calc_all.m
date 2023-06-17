%% Decoding analyses for stimuli in all durations (300,900,1500ms):
% Result section 4 (Figure 4 & Supp Figure 7)
% Take the output from here to decoding_fig_all.m

clear all;clc;close all

%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));

%% settings

% Define what segments we load (from specific_segs.m - run that before):
% Subjects:
subs        = 4:10; % 1:3 or 4:10 (need to run both)

% Categories: (face,watch,object,animal)
cat_nums    = 1:2;  % for almost everything, except Figure S7d which requires [1,3]

% Define what electrodes\ROIs we use to do the decoding:
%   (all non-noisy electrodes & all ROIs are loaded in the initial segments)

% which electrodes? 'resp' (responsive)\'pos' (responsive & increasing from baseline) [can also use 'all', but it wasn't in the paper]
%   Use 'resp' for almost everything
%   'pos' for Figure S7e
e_subsample = 'resp';

% ROIs: (Occ,VT,Par,PFC)
ROIs        = 1:4; % just run this for both patient groups
% in the main we use regions 1,2,4 from subj 4:10 and region 3 from 1:3,
% supp also includes regions 1:2 from subj 1:3

%% load data (created with specific_segs.m):

load([DATA_FOLDER,'figure_settings.mat'],'ROInames')
% constant settings
d_str = 'all'; downsample_hz = 200; n_rep = 1; av_in_cat = false;
segs_desc = sprintf('%drep_dur%s_subs%s_%dHz_cat%s_av%s',n_rep,d_str,join(string(subs),''),...
    downsample_hz,join(string(cat_nums),''),string(av_in_cat));
load([DATA_FOLDER, sprintf('segs_%s.mat',segs_desc)]);

% some more followups on the settings
tmp = tabulate(categories_vec); if any(tmp(:,2)<10); n_folds = 4; else; n_folds = 5; end
switch e_subsample
    case 'resp'
        e_filt = elec_info.IsResponsive;
    case 'pos'
        e_filt = elec_info.IsResponsive & elec_info.RespSign==1;
end
segs = segs(:,e_filt,:); elec_info = elec_info(e_filt,:);

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
cfg_decoding.k          = n_folds;

% statistics settings: see decoding_stats.m for more details
cfg_stats = [];
% stat_type: 'analytic'\'perm'\'cluster'\'max', how to do correction for multiple comparisons
%   cluster uses cluster-based permutations, max with max-statistic control for MC
%   the other two are corrected with Benjamini-Hochberg FDR
%   note I added the max option only after completing the manuscript so it wasn't tested thoroughly
cfg_stats.stat_type  = 'cluster'; % I use cluster here and compute FDR when needed in the figures code for matrices etc.
cfg_stats.n_sides    = 1; % !!! right sided (larger than)\dual sided test
cfg_stats.p_thresh   = 0.05; % threshold for significance (after fdr\for the clusters)
cfg_stats.n_perm     = 1000;
cfg_stats.acc        = strcmp(cfg_decoding.metric,'acc');
% cluster specific settings (see inside the function)
cfg_stats.n_clusters = Inf;
cfg_stats.firstlevel_type   = 'stat';
cfg_stats.firstlevel_thresh = 0.1; % note that the thresh is relative to chance level!

% the split into different durations is done in a dedicated function, this
% is the input structure for that:
add_inf = [];
add_inf.cfg_decoding   = cfg_decoding;
add_inf.cfg_stats      = cfg_stats;
add_inf.mask_data      = true; % tells stats to actually run

% we repeat the statistics with point by point permutations here which was
% used for some of the figures and then we don't need to save the permutations themsevles
cfg_stats_perm           = cfg_stats;
cfg_stats_perm.stat_type = 'perm';

%% calculation (can take some time)
% most of the action is in dec_mult_dur.m
% For each comparison we compute the temporal generalization matrix (King &
%   Dehaene, 2014), and the 'standard' decoding is the diagonal of this.
%   stats are computed for both so we won't 'over-correct' for multiple
%   comparisons if we don't use the full temporal generalization matrix.

for reg_i = ROIs
if reg_i<=6; ROI_idx = grp2idx(elec_info.ROI); else; ROI_idx = grp2idx(elec_info.ROISplit); end % use ROIsplit if any of regions 7-10 requested
segs_roi = segs(:,ROI_idx==reg_i,:); n_elec = size(segs_roi,2);
if n_elec < 1; continue; end % nothing to run on (e.g. 1:3 PFC there are no responsive electrodes)

rng(10) % for reproducability
% some fixes for how I presented the data in the paper
data_out = dec_mult_dur(segs_roi, time_vec, {'categories', categories_vec, cat_names}, add_inf);

% add some more stats to the save and remove permutations to make it less heavy
% rerunning just the diagonals with cluster perm & point by point + FDR,
% and both only using the time between the two durations as done in the
% manuscript (this is not what I used for the matrices, so that stays the same)
res_diag_masks_perm = nan(size(data_out.res_diag_masks)); % point by point + FDR 
diff_diag_masks_perm = nan(size(data_out.diff_diag_masks));
n_dur = length(data_out.dur_unique);
diff_diag_masks_clust2 = nan(size(data_out.diff_diag_masks)); pvals_clust2 = cell(1,n_dur-1);
for dur = 1:n_dur
    keep_idx = time_vec <= data_out.dur_unique(dur)+600;
    vals_cur = diag(data_out.res_mats(keep_idx,keep_idx,dur)); perms_cur = get_diags(data_out.perms(keep_idx,keep_idx,:,dur));
    res_diag_masks_perm(keep_idx,:,dur) = decoding_stats(cfg_stats_perm, vals_cur, perms_cur);
    if dur < 3
        keep_idx = keep_idx & time_vec >= data_out.dur_unique(dur);
        vals_cur = diag(diff(data_out.res_mats(keep_idx,keep_idx,dur:(dur+1)),[],3));
        perms_cur = get_diags(diff(data_out.perms(keep_idx,keep_idx,:,dur:(dur+1)),[],4));
        diff_diag_masks_perm(keep_idx,:,dur) = decoding_stats(cfg_stats_perm, vals_cur, perms_cur);
        [diff_diag_masks_clust2(keep_idx,:,dur), pvals_clust2{dur}]= decoding_stats(data_out.cfg_stats_diff, vals_cur, perms_cur); % redoing because of the change in the time-points
    end
end
fields = {'res_diag_masks_perm','diff_diag_masks_perm','diff_diag_masks_clust2','pvals_clust2'};
for f = 1:length(fields); eval(sprintf('data_out.%s = %s;', fields{f}, fields{f})); end
data_out.res_mats = 100*data_out.res_mats; % turn it to % success (not decimal out of 1)
data_out.perms    = 100*data_out.perms;
data_out = rmfield(data_out,'perms'); % can remove this part to keep the permutations
save([DATA_FOLDER,sprintf('decoding_%s_%s.mat',ROInames(reg_i),general_save_name)],'data_out') % can also add segs_roi for completeness and keep perms etc.
end

% saving settings relevant for all ROIs (to be used later with the figures)
add_inf.cfg_stats.firstlevel_thresh = 10; % adjust to this units change
data_out.cfg_stats_diff.firstlevel_thresh = 10; % adjust to this units change
tmp = split(cat_names',' '); categ_only = unique(tmp(:,1));
comp_name      = join(cellfun_wrap(@(x) [upper(x(1)),x(2:end)],categ_only),'-'); comp_name = comp_name{:};
save([DATA_FOLDER, sprintf('decoding_settings_%s.mat', general_save_name)],'time_vec','add_inf','comp_name');
% you can add categories_vec, just not used in the figure generating code.