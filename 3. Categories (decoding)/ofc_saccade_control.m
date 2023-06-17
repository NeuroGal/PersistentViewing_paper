%% Comparison of decoding in OFC with\without saccades (Figure S6e)
% Take the output from here to decoding_fig_long.m

clear all;clc;close all
%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));

%% Start by loading the relevant data:
% 'wins'           - the time-windows we are working with (4x2, each row is [window_start,window_end])
% 'segs_wins'      - S8 responses averaged in the four bins (trials x electrodes x 4-time-windows).
%                    Electrodes already includes only reponsive OFC electrodes on grids fully in OFC (as described in the methods).
%                    Trials include only images in the 4 analyzed categories over 900 ms.
% 'cat_stim_mat'   - matrix of n_trials x 2 with [category, stimulus_id] in each row.
% 'saccade_in_win' - n_trials x 4-time-windows with True\False in each position indicating whether
%                    a saccade was detected in that window or not.%

load([DATA_FOLDER,'segs_ofc_saccade_control.mat'])

% 'saccade_in_win' was computed using the following code:
%   E_data = ... % time x 1, LFP data from the relevant temporal pole electrode in S8 which contains saccadic spike potentials (STAR methods) 
%                % this is not currently on OSF, but please contact us if you are interested in running this analysis.
%   C = filtSRP(double(E_data),SDATA.info.sampRate); % run match template filter from Keren, Yuval-Greenberg and Deouell, 2010
%                                                    % (https://doi.org/10.1016/j.neuroimage.2009.10.057) - also included in the code folder
%   sac = [0 ; diff(C>3*std(C))>0];
%   sac_seg = ... % trials x time, segmented 'sac' around stimulus onset excluding time of seizure activity as done for the main segments
%                 % and including only non-target stimuli in the relevant categories over 900 ms
%   saccade_in_win = cellfun_wrap(@(x) sum(sac_seg(:,ismember(time_vec, x(1):x(2))),2), num2cell(wins,2)',true)>0;

%% Settings

cfg_decoding = [];
cfg_decoding.preprocess = 'undersample';
cfg_decoding.metric = 'auc'; 

n_perm = 1000; 
win_names = join(string(wins),'-')+"ms";
cat_names = {'Face','Watch','Object','Animal'};
n_comp = nchoosek(4,2); % the comparison order and names are the same as they are in the standard case, so not saving here
n_win = length(win_names);

% this is only used so decoding_wrapper will export permutations, since we use 2 windows (for convenience, see below) 
% when we don't want that to be FDR corrected we do the stats outside (Bonferroni correction is done in the figure code).
cfg_stats = []; cfg_stats.stat_type = 'perm'; cfg_stats.n_sides = 1;
cfg_stats.n_perm = n_perm; cfg_stats.p_thresh = 0.05; cfg_stats.acc=false; cfg_stats.chance_level = 0.5;

%% Find how many saccades occured for each image (II in the Figure)

[~,unique_idx,stim_idx]=unique(cat_stim_mat,'rows');
unique_image_sac = nan(length(unique_idx), size(wins,1));
for idx = 1:length(unique_idx)
    unique_image_sac(idx,:) = 1*all(saccade_in_win(stim_idx==idx,:),1)+1*any(saccade_in_win(stim_idx==idx,:),1); % 0-no sac,1-both no and yes sac, 2-sac all trials
end
cat_nums_tab = tabulate(cat_stim_mat(unique_idx,1));
tmp = cellfun_wrap(@(x) tabulate(x), mat2cell(unique_image_sac,cat_nums_tab(:,2),[1 1 1 1])); % category x window
trial_counts = permute(cellfun_wrap(@(x) permute(x(:,2),[2 3 1]), tmp, true),[3 2 1]); % type (no sac, part w sac part w\o, all sac) x window x category

%% run the decoding calculation (I,III,IV in the Figure)

perms_all = nan(n_comp,n_win,2,n_perm); perf_all = nan(n_comp,n_win,2); % 2 types - all trials, no sac
for type_idx = 1:2
for w_idx = 1:n_win
if type_idx == 1
    trial_filt = true(size(saccade_in_win(:,w_idx)));
elseif type_idx == 2
    trial_filt = ~saccade_in_win(:,w_idx);
end
rel_segs_wins = segs_wins(trial_filt, :, :);
[cat_stim_unique,~,unique_idx] = unique(cat_stim_mat(trial_filt,:), 'rows');
segs_wins_unique = nan([size(cat_stim_unique,1), size(rel_segs_wins,[2 3])]);
for stim_i = 1:size(cat_stim_unique,1) % merge repetitions of the same stimulus
    segs_wins_unique(stim_i, :, :) = mean(rel_segs_wins(unique_idx==stim_i,:,:),1);
end
decoded_segs = segs_wins_unique(:,:,[w_idx 1]); % adding the first window in all cases just so the decoding_wrapper will work without issues (to have something in the time dimension)
categories_vec = cat_stim_unique(:,1);

% as explained in the main decoding_calc_long.m file - originally the
% categories were ordered differently, to get the same numbers I am
% flipping it and will flip again below:
categories_vec_tmp = grp2idx(categorical(categories_vec,[1,3,4,2]));
[categories_vec_tmp,tmp_idx]=sort(categories_vec_tmp);
decoded_segs_tmp = decoded_segs(tmp_idx,:,:); reorder_idx = [3 1 2 5 6 4];

rng(10) % for reproducability
[perf, ~, ~, perms] = decoding_wrapper(decoded_segs_tmp, categories_vec_tmp, cfg_decoding, cfg_stats);
perf_all(:,w_idx,type_idx) = 100*squeeze(perf(1,1,reorder_idx)); % extract only the relevant window
perms_all(:,w_idx,type_idx,:) = squeeze(100*cellfun_wrap(@(x) x(1,1,:), perms(reorder_idx), true)); % same
end
end
pvals_all = (sum(perf_all<=perms_all,4)+1)/(n_perm+1);
save([DATA_FOLDER,'decoding_ofc_saccade_control.mat'],'perf_all','wins','win_names','cat_names','pvals_all','trial_counts')