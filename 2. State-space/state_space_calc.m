%% Analyses for results section 2 (Figure 1 & Supp Figure 2-3)
% take the output from here to state_space_fig.m - run once with
% mult_dur_mode = true and once = false.
% notice this uses max_stat_correction from the time_resolves_stats toolbox

clear all;clc;close all
%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));
load([DATA_FOLDER,'figure_settings.mat'],'ROInames_full')

% settings:
mult_dur_mode = false; % False: Figure 2a, Supp Figure 2, Supp Figure 3e; True: Figure 2b, Supp Figure 3a-d
subs   = 1:10;
if ~mult_dur_mode
d_str = 'long';  cat_nums  = 1:4; stim_durs = 900; dur_pairs = [1 1];
else
d_str = 'all'; cat_nums  = 1:2; stim_durs = 300:600:1500; dur_pairs = [1 2; 2 3];
% comparing only 1-2 and 2-3 to save time, we can also add 1-3 (or just use: nchoosek(1:n_dur, 2))
end
n_dur  = length(stim_durs); n_pair = size(dur_pairs,1); pair_times = stim_durs(dur_pairs);
n_perm = 1000;
sample_hz = 1000;
if strcmp(d_str,'long');tv_end = 900;elseif strcmp(d_str,'all'); tv_end = 2100; end
stat_info = []; stat_info.p_thresh = 0.05; stat_info.n_sides = 1; stat_info.chance_level = 0;

% usefull functions for later
normer = @(dat) squeeze(cellfun(@norm, num2cell(dat,2))); % assumes elec is dim 2
normalizer = @(main_array, perm_array, perm_dim) (main_array-mean(perm_array,perm_dim))./std(perm_array,[],perm_dim);

%% Load basic segs & get segs from each subj - takes a few min to load, but then the calculation is very fast

if ~exist('segs_basic','var'); load([DATA_FOLDER, sprintf('basic_segs_%dHz.mat',sample_hz)]);cat_names_old = cat_names; elec_info_old = elec_info; end
t_inds   = time_vec_basic > -100 & time_vec_basic <= tv_end; % define the output time range
time_vec = time_vec_basic(t_inds); n_time = length(time_vec);

segs_per_sub = cell(length(subs),1); info_per_sub = cell(length(subs),1);
for s = subs
    trial_info = trial_info_basic{s};
    % select relevant categories
    include = ismember(trial_info.cat_nums, cat_nums);
    % select relevant durations
    if ~mult_dur_mode
        cat_stim_mat = [trial_info.cat_nums, trial_info.stim_ids];
        varnames = {'categ','stim_id'};
        include = include & trial_info.durations >= tv_end;
    elseif mult_dur_mode
        cat_stim_mat = [trial_info.cat_nums, trial_info.stim_ids, trial_info.durations];
        varnames = {'categ','stim_id','dur'};
        include = include & ismember(trial_info.durations,[300 900 1500]);
    end
    % select relevant electrodes
    e_filt = elec_info_old.IsResponsive(elec_info_old.Patient==s);
    segs = double(segs_basic{s}(include,e_filt,t_inds));
    cat_stim_mat = cat_stim_mat(include, :);
    
    % average repetitions of the same stimulus
    [cat_stim_unique_mat,~,unique_idx] = unique(cat_stim_mat, 'rows');
    stim_info = array2table(cat_stim_unique_mat,'VariableNames',varnames);
    segs_choice = nan([size(stim_info,1), size(segs,[2 3])]);
    for stim_i = 1:size(stim_info,1)
        segs_choice(stim_i, :, :) = mean(segs(unique_idx==stim_i,:,:),1);
    end
    segs = segs_choice;
    if mult_dur_mode % insert nan where another stimulus may have starter (600 ms is the shortest ISI)
        for dur = stim_durs; segs(stim_info.dur == dur, :, time_vec > dur+600) = nan; end
    end
    segs_per_sub{subs==s} = segs;
    info_per_sub{subs==s} = [stim_info,array2table(repmat(s,height(stim_info),1),'VariableNames',{'patient'});];
    fprintf('Done with subj #%d\n',s)
end
cat_names = cat_names_old(cat_nums); elec_info = elec_info(elec_info.IsResponsive,:); % final fixes
stim_info_all = cat(1,info_per_sub{:}); % get the stim info to one table for all subjects
elec_per_subj = cellfun(@(x) size(x,2), segs_per_sub); % number of electrodes per subjects
n_elec_total  = sum(elec_per_subj);

% get unique category (similar to category_vec elsewhere)
if ~mult_dur_mode
    stim_info_all.unique_cat = grp2idx(categorical(stim_info_all.categ));
elseif mult_dur_mode
    stim_info_all.unique_cat = grp2idx(categorical(stim_info_all.categ*100 + stim_info_all.dur/100));
    cat_names = cellstr(string(cat_names) + ' ' + string(stim_durs)');  cat_names = cat_names(:); % order cat1 300, cat1 900, cat1 1500 then cat2 etc.
end 
n_cat = length(unique(stim_info_all.unique_cat));
ROI_idx = grp2idx(elec_info.ROI); n_reg = 6;

%% CI to the means per category by jackniffing + follow up calculations ~10-15 min for one duration, double for mult dur

run_ci = true;
% currently doesn't run speed calculations for mult dur (not used in the paper) - you can change that
if run_ci
%% segments used to calculate distance to baseline, speed of the trajectory etc
% in particular: confidence intervals for supp figure 2b,d
rng(10) % for reproducability
segs_ci = nan(n_cat, n_elec_total, n_time, n_perm + 1);
elec_ind = 0;
for s=subs
    segs = segs_per_sub{subs==s}; n_elecs = size(segs,2);
    sub_unique_cat = stim_info_all.unique_cat(stim_info_all.patient == s);
    for perm_p = 1:(n_perm+1)
        for categ = 1:n_cat
            ids = find(sub_unique_cat==categ);
            if perm_p > 1; ids = randsample(ids,length(ids)-1); end  %jackknife
            segs_ci(categ, (elec_ind+1):(elec_ind+n_elecs), :, perm_p) = mean(segs(ids,:,:),1);
        end
        if mod(perm_p,100)==1; fprintf('Done perm %d/%d\n',perm_p-1,n_perm); end
    end
    elec_ind = elec_ind + n_elecs;
    fprintf('Done with subj #%d\n',s)
end

%% folloup calculations on that (baseline_dist, speed etc) (a few min)
% Calculating baseline_dist
%   Used for supp 2b (and plotted also in all of the trajectory plots as an inset, but not with CI)
%   assumes the segments were baselined, otherwise 'normer' needs to take segs-baseline as input
dims = [n_cat, n_time, n_reg, n_perm + 1];
baseline_dist = nan(dims);
for reg_i = 1:n_reg
    baseline_dist(:,:,reg_i,:) = normer(segs_ci(:, ROI_idx==reg_i, :,:)); 
    fprintf('Done with %s (reg #%d)\n',ROInames_full(reg_i),reg_i)
end

% Calculating speed 
%   Used for Figure S2c-d
%   the distance between the curr tp and the one before (technically should
%   be normalized by the time between 2 tp, but it's 1 second for us)
if ~mult_dur_mode
speed = nan(dims);
for reg_i = 1:n_reg
    speed(:,2:end,reg_i,:) = normer(diff(segs_ci(:, ROI_idx==reg_i, :,:),[],3)); 
    fprintf('Done with %s (reg #%d)\n',ROInames_full(reg_i),reg_i)
end

% baseline_dist and speed CI peaks and attenuation (pretty fast)
% used for Figure S2b,d
t_zero = find(time_vec == 0);
dims = [n_cat, n_reg, n_perm + 1, 2]; % last is baseline_dist, speed
peak_time = nan(dims); percent_atten = nan(dims); 
for meas = 1:2
    if meas == 1; rel_array = baseline_dist; elseif meas == 2; rel_array = speed; end
    rel_array = permute(rel_array(:,(t_zero+1):end,:,:), [1 3 2 4]); % now it's cat,roi,time,perm and only after t_zero
    [peak_val, peak_idx] = max(rel_array, [], 3);
    peak_time(:, :, :, meas) = time_vec(peak_idx + t_zero);

    dat_end = mean(rel_array(:, :, time_vec((t_zero+1):end)>=800 & time_vec((t_zero+1):end)<=900, :),3);
    percent_atten(:, :, :, meas) = 100*(peak_val-dat_end)./peak_val; % not comparing the actual values because they depend on the number of electrodes and more
    fprintf('Done with meas %d/2\n', meas)
end
% Speed around peak of baseline_dist (pretty fast)
% Used in Figure S2d
t_rng = 100; % in ms (though with 1000Hz it's the same)
speed_finder = @(vec ,t0) (mean(vec(t0>time_vec & time_vec>=t0-t_rng))-mean(vec(t0<time_vec & time_vec<=t0+t_rng)))/vec(time_vec == t0); % doesn't have to be speeds...
speed_diff = cellfun(@(vec, t0) 100*speed_finder(vec,t0), squeeze(num2cell(speed,2)), num2cell(peak_time(:,:,:,1)));
end
end

%% Permutation of category affiliation + category selectivity index - <10min on my computer (not for mult dur)
% calculating multivariate category selectivity + stats (supp figure 3e)

if ~mult_dur_mode
rng(10) % for reproducability
segs_shuff = nan(n_cat, n_elec_total, n_time, n_perm + 1); 
elec_ind = 0;
for s=subs
    segs = segs_per_sub{subs==s}; n_elecs = size(segs,2);
    sub_unique_cat = stim_info_all.unique_cat(stim_info_all.patient == s);
    for perm_p = 1:(n_perm+1)
        ids = sub_unique_cat;
        if perm_p > 1; ids = ids(randperm(length(ids))); end
        for categ = 1:n_cat
            segs_shuff(categ, (elec_ind+1):(elec_ind+n_elecs), :, perm_p) = mean(segs(ids==categ,:,:),1);
        end
        if mod(perm_p,100)==1; fprintf('Done perm %d/%d\n',perm_p-1,n_perm); end
    end
    elec_ind = elec_ind + n_elecs;
    fprintf('Done with subj #%d\n',s)
end
% Calculating cat variability - based on segs_shuff
calc_multivar_catsel = @(array) sqrt(mean(normer(array-mean(array,1)).^2,1));
multivar_catsel_shuff = nan(n_time, n_reg, n_perm + 1);
for reg_i = 1:n_reg
    multivar_catsel_shuff(:, reg_i, :) = calc_multivar_catsel(segs_shuff(:, ROI_idx==reg_i, :, :));
    fprintf('Done with %s (reg #%d)\n',ROInames_full(reg_i),reg_i)
end
% normalizing by the permutation distribution!
multivar_catsel_shuff = normalizer(multivar_catsel_shuff, multivar_catsel_shuff(:,:,2:end), 3);

max_stat_array = permute(multivar_catsel_shuff,[1 4 3 2]); % move the dimensions we don't want to correct to the end so it will work with the format of max_stat_correction.
[catsel_masks, ~, catsel_thresh] = max_stat_correction(max_stat_array(:,:,1,:), max_stat_array(:,:,2:end,:), stat_info);
catsel_masks = squeeze(catsel_masks); multivar_catsel = multivar_catsel_shuff(:,:,1); % save only what we need for plotting
end
 
%% Distance between durations & permutation of duration affiliation - 30-40 min (mult dur longer segments)

if mult_dur_mode
%% CI bet-durs distances < 10 min
% Used in Supp Figure 3c-d
%   distances from other stim of your category, at the same time
%   only until 1500 since we use differences later anyway (the rest will be nan)
%   notice this is normalized by the null distribution below!
dist_bet_durs_ci = nan([(n_cat/n_dur)*n_pair, sum(time_vec<=1500), n_reg, n_perm + 1]); 
long_ids  = dur_pairs(:,2) + (0:n_dur:(n_cat-n_dur)); long_ids = long_ids(:);
short_ids = dur_pairs(:,1) + (0:n_dur:(n_cat-n_dur)); short_ids = short_ids(:);
subt_segs = segs_ci(short_ids,:,time_vec<=1500,:)-segs_ci(long_ids,:,time_vec<=1500,:);
for reg_i = 1:n_reg
    dist_bet_durs_ci(:,:,reg_i,:) = normer(subt_segs(:,ROI_idx==reg_i,:,:)); 
    fprintf('Done with %s (reg #%d)\n',ROInames_full(reg_i),reg_i)
end
clear subt_segs

%% Permutation of duration affiliation <10min
%   Shuffling durations and calculating the mean of the shuffle,
%   used below to calculate the permutation null dist for distance between
%   them and difference of the dist-to-baseline

rng(10) % for reproducability [originally I ran also the 300-1500 comparison so it's going to be very similar but not identical]
dur_shuff_sub = nan([(n_cat/n_dur)*n_pair, n_elec_total, sum(time_vec<=1500), n_perm + 1, 2]);
% dim1 - for each category (n_cat/n_dur) we loop over n_pairs, last-dim - dur1, dur2 in pair
cat_nums = unique(stim_info_all.categ); % extra precaution if the settings in the beginning change (e.g. use categories 1,3)
elec_ind = 0;
for s=subs
    segs = segs_per_sub{subs==s}(:, :, time_vec<=1500); n_elecs = size(segs,2);
    for categ = 1:(n_cat/n_dur)
        cat_filt = info_per_sub{subs==s}.categ == cat_nums(categ); % the actual category (not taking each dur separately)
        for pair = 1:n_pair
            dur_filt = ismember(info_per_sub{subs==s}.dur,stim_durs(dur_pairs(pair,:)));
            cur_segs = segs(cat_filt & dur_filt, :, :);
            sub_dur = info_per_sub{subs==s}.dur(cat_filt & dur_filt);
            for perm_p = 1:(n_perm+1)
                ids = sub_dur;
                if perm_p > 1; ids = ids(randperm(length(ids))); end
                for dur_idx = 1:2
                    dur_shuff_sub(n_pair*(categ-1)+pair, (elec_ind+1):(elec_ind+n_elecs), :, perm_p, dur_idx) = mean(cur_segs(ids==stim_durs(dur_pairs(pair,dur_idx)),:,:),1);
                end
                if mod(perm_p,100)==1; fprintf('Done perm %d/%d\n',perm_p-1,n_perm); end
            end
        end
    end
    elec_ind = elec_ind + n_elecs;
    fprintf('Done with subj #%d\n',s)
end
clear segs cur_segs

%% Compute the shuffled null distributions ~15-20min
%   Used in Figure 2b and Figure S3a-d for the statistics
%   Shuffled null dist for distance between trajectories of different durations\ difference in dist-to-baseline

dims = [(n_cat/n_dur)*n_pair, sum(time_vec<=1500), n_reg, n_perm + 1];
dist_bet_durs_shuff = nan(dims); baseline_dist_diff_shuff = nan(dims);
for reg_i = 1:n_reg
    segs_use = dur_shuff_sub(:,ROI_idx==reg_i,:,:,:);
    dist_bet_durs_shuff(:,:,reg_i,:) = normer(diff(segs_use,[],5)); 
    fprintf('Dist between durs shuff - done with %s (reg #%d)\n',ROInames_full(reg_i),reg_i)
    baseline_dist_diff_shuff(:,:,reg_i,:) = normer(segs_use(:,:,:,:,2))-normer(segs_use(:,:,:,:,1)); 
    fprintf('Baseline dist diff shuff - done with %s (reg #%d)\n',ROInames_full(reg_i),reg_i)
end
clear segs_use dur_shuff_sub

%% normalizing by the permutation distribution, stats & fixing things to the correct dimensions for plotting
dist_bet_durs_ci = normalizer(dist_bet_durs_ci, dist_bet_durs_shuff(:,:,:,2:end), 4);
dist_bet_durs_shuff = normalizer(dist_bet_durs_shuff, dist_bet_durs_shuff(:,:,:,2:end), 4);
baseline_dist_diff_shuff = normalizer(baseline_dist_diff_shuff, baseline_dist_diff_shuff(:,:,:,2:end), 4);

tmp = nan(size(baseline_dist)); tmp(short_ids,time_vec<=1500,:,:) = dist_bet_durs_ci; dist_bet_durs_ci = tmp;
tmp = nan(size(baseline_dist)); tmp(short_ids,time_vec<=1500,:,:) = dist_bet_durs_shuff; dist_bet_durs_shuff = tmp;
tmp = nan(size(baseline_dist)); tmp(short_ids,time_vec<=1500,:,:) = baseline_dist_diff_shuff; baseline_dist_diff_shuff = tmp;

dist_bet_durs_masks = false(size(dist_bet_durs_shuff,1:3)); 
baseline_dist_diff_masks = false(size(baseline_dist_diff_shuff,1:3)); 
for pair = 1:n_pair % different time points
rel_comps = dur_pairs(pair,1) + n_dur*(0:(n_cat/n_dur - 1));
rel_time = pair_times(pair, 1) < time_vec & time_vec <= pair_times(pair, 2); 
max_stat_array = permute(dist_bet_durs_shuff(rel_comps,rel_time,:,:),[2 5 4 1 3]); % move the dimensions we don't want to correct to the end so it will work with the format of max_stat_correction.
dist_bet_durs_masks(rel_comps,rel_time,:) = permute(max_stat_correction(max_stat_array(:,:,1,:,:), max_stat_array(:,:,2:end,:,:), stat_info),[4 1 5 2 3]);
max_stat_array = permute(baseline_dist_diff_shuff(rel_comps,rel_time,:,:),[2 5 4 1 3]); % same
baseline_dist_diff_masks(rel_comps,rel_time,:) = permute(max_stat_correction(max_stat_array(:,:,1,:,:), max_stat_array(:,:,2:end,:,:), stat_info),[4 1 5 2 3]);
end
dist_bet_durs = dist_bet_durs_ci(:,:,:,1); % only this is needed for the plot
baseline_dist = baseline_dist(:,:,:,1); % same
end

%% SAVE

save_in_both = {'time_vec','n_reg','cat_names','baseline_dist'};
if ~mult_dur_mode
    mean_segs = segs_ci(:,:,:,1);
    save([DATA_FOLDER,'state_space_calc.mat'],save_in_both{:},'mean_segs','ROI_idx','speed','peak_time','percent_atten','speed_diff','catsel_masks','catsel_thresh','multivar_catsel');
else
    save([DATA_FOLDER,'state_space_calc_mult_dur.mat'],save_in_both{:},'stim_durs','dur_pairs','pair_times','baseline_dist_diff_masks', 'dist_bet_durs', 'dist_bet_durs_masks')
end