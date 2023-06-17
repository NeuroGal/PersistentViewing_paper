%% Get from the basic segments to specific segments needed for each analysis

clear all;clc;close all
% DATA_FOLDER = pwd; % change to the folder where you put basic_segs_XHz.mat and the other data
% CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
% addpath(genpath(CODE_FOLDER));

%% settings (see more details in each of the analysis scripts)

% Sampling rate - use the 200 Hz or 1000 Hz version of basic_segs.m?
%   Use 200 Hz for decoding and single-exemplar analyses, 1000 Hz for all else.
sample_hz       = 1000;

% Which patients to use?
%   Main Figures - Figures 1-3: subj 1:10, Figure 4: subj 4:10 and 1:3, Figure 5: subj 1,2,5,7,8
%   Supplementary uses these settings or single patients (Figures S5-6) \ patients [8,10] (Figure S6e only) 
%       and [1:2,5:10] for some of the single-exemplar controls
subs            = 1:10;

% Which durations to use?
%   'long' - merge all stimuli >= 900ms, and then the time-vector is cropped at 900ms. 
%   'all'  - keeps all stimuli durations distinct.
%   Use 'long' for everything except Figure 4 & Supp Figure 7
d_str           = 'long';   

% How many repetitions?
%   1 - images which were viewed by all patients at least once are kept (if
%       tv_end = 900 this doesn't have to be the same duration, only that all
%       patients viewed them for durations >= 900ms. tv_end = 2100 ensures they
%       were all viewed in the same duration.
%   2 - images viewed at least twice.
% Use 2 for single-exemplar analyses, 1 for all else
n_rep           = 1;

% Which categories to keep? (order: face,watch,object,animal)
%   Use 1:4 for most things, [1,2] or [1,3] for specific decoding analyses
%   (explained in the _calc.m scripts)
cat_nums        = 1:4;

% Average all stimuli together by categories?
%   False: keeps different stim_id separate
%   True:  does the merge
% Use for Figures 1-2, supp Figures 1-3, otherwise false
av_in_cat       = true; 

% Output:
% 'segs'      - stimulus x electrode x time, similar to segs_basic but merged across patients
% 'stim_info' - table with fields:
%       'categ'   - category, 1:n_cat (so even if you chose [1,3] above this will be 1:2, the categories are named in cat_names.
%       'stim_id' - unique stimulus id, changed to 999 if av_in_cat
%       'dur'     - duration of stimulus, only added if tv_end=2100 (otherwise it's all >=900ms)
%       'rep'     - repetition (1/2), only when n_rep = 2.
% 'time_vec'
% 'categories_vec' - identical to stim_info.categ when it's only the long
%                    durations merged, but when all durations are kept
%                    separately this is a unique identified merging
%                    categ & dur.
% 'cat_names' - cell array of category titles, includes the duration in the
%               case of all stimuli (corresponding to categories_vec).
% 'elec_info' - similar to the structure in segs_basic.m, but only with the
%               relevant electrodes (corresponding to the chosen subjects).

%% load data from basic_segs_XXHz.mat

load([DATA_FOLDER, sprintf('basic_segs_%dHz.mat',sample_hz)]); % notice elec_info & cat_names are changed below so you need to load them again to run with other settings
% segs_basic excludes already trials including ictal spikes\ other HF noise
% after smoothing by 50ms moving average and baseline subtraction of -300 to 0 ms pre-stim.

segs_desc = sprintf('%drep_dur%s_subs%s_%dHz_cat%s_av%s',n_rep,d_str,join(string(subs),''),sample_hz,join(string(cat_nums),''),string(av_in_cat))
if strcmp(d_str,'long');tv_end = 900;elseif strcmp(d_str,'all'); tv_end = 2100; end

% get elec_info for the electrodes we selected (based on the subjects)
elec_info = elec_info(ismember(elec_info.Patient,subs),:);
% define the relevant time range and categories
t_inds   = time_vec_basic > -100 & time_vec_basic <= tv_end; % define the output time range
time_vec = time_vec_basic(t_inds); cat_names = cat_names(cat_nums);

%% Get to specific segs in each subject

segs_specific = cell(length(subs),1); info_specific = cell(length(subs),1);
for s = subs
    segs = segs_basic{s}; trial_info = trial_info_basic{s};

    % select relevant categories
    include = ismember(trial_info.cat_nums, cat_nums);
    
    % select relevant durations
    switch d_str
        case 'long'
            cat_stim_mat = [trial_info.cat_nums, trial_info.stim_ids];
            varnames = {'categ','stim_id','n_trials'};
            include = include & trial_info.durations >= tv_end;
        case 'all'
            cat_stim_mat = [trial_info.cat_nums, trial_info.stim_ids, trial_info.durations];
            varnames = {'categ','stim_id','dur','n_trials'};
            include = include & ismember(trial_info.durations,[300 900 1500]);
    end
    segs = segs(include, :, :); cat_stim_mat = cat_stim_mat(include, :);
    
    % average repetitions of the same stimulus (to 1/2 repetitions per image)
    [cat_stim_unique_mat,~,unique_idx] = unique(cat_stim_mat, 'rows');
    tmp = tabulate(unique_idx); stim_ids_n = tmp(:,2);
    stim_info = array2table([cat_stim_unique_mat, stim_ids_n],'VariableNames',varnames);
    if n_rep==1
        segs_choice = nan([size(stim_info,1), size(segs,[2 3])]);
        for stim_i = 1:size(stim_info,1)
            segs_choice(stim_i, :, :) = mean(segs(unique_idx==stim_i,:,:),1);
        end
    elseif n_rep==2
        enough_reps = find(stim_ids_n >= n_rep); stim_info = stim_info(enough_reps, :);
        segs_choice = nan([size(stim_info,1), size(segs,[2 3]), 2]);
        for stim_i = 1:size(stim_info,1)
            rel_idx = find(unique_idx==enough_reps(stim_i));rel_idx = rel_idx(randperm(length(rel_idx)));
            split_idx = floor(stim_info.n_trials(stim_i)/2);
            part1 = mean(segs(rel_idx(1:split_idx),:,:),1);
            part2 = mean(segs(rel_idx((split_idx+1):end),:,:),1);
            segs_choice(stim_i, :, :, :) = cat(4, part1, part2);
        end
    end
    segs = segs_choice(:,:,t_inds,:); % crop the segments to the relevant time indices
    
    if av_in_cat % if we want to average all images of each category per electrode
        rel_info = stim_info(:,ismember(stim_info.Properties.VariableNames,{'categ','dur'}));
        [~,idx_to_unique,idx_from_unique] = unique(rel_info,'rows','stable');
        dims = size(segs); dims(1) = length(unique(idx_from_unique));
        segs_av = nan(dims);
        for i = 1:dims(1)
            segs_av(i, :, :, :) = mean(segs(idx_from_unique==i,:,:,:),1);
        end
        stim_info = stim_info(idx_to_unique, :);
        tmp = tabulate(idx_from_unique); stim_info.n_trials = tmp(:,2); 
        stim_info.stim_id = repmat(999,length(stim_info.stim_id),1);
        segs = segs_av;
    end
    segs_specific{subs==s} = segs; info_specific{subs==s} = [stim_info,array2table(repmat(s,height(stim_info),1),'VariableNames',{'patient'});];
    fprintf('Done with subj #%d\n',s)
end

%% merge across subjects

info_specific_all = cat(1,info_specific{:}); % get the stim info to one table for all subjects
elec_per_subj = cellfun(@(x) size(x,2), segs_specific); % number of electrodes per subjects
% find the images which all subjects have seen
cat_stim_ids_idx = ismember(info_specific_all.Properties.VariableNames,{'categ','stim_id','dur'});
[cat_id_unique, ~, to_unique] = unique(info_specific_all(:,cat_stim_ids_idx), 'rows');
tmp = tabulate(to_unique); value_count = tmp(:,2);
stim_info = cat_id_unique(value_count >= length(segs_specific), :);

% merge segments from all subjects
seg_dims = [length(stim_info.categ),sum(elec_per_subj),length(time_vec), n_rep]; segs = nan(seg_dims);
elec_ind = 0;
n_trials_per_sub = nan(size(stim_info,1),length(subs)); % needed later in some cases (e.g. to create the overall mean in univariate calc)
for s=subs
    [~,ids_in_global,ids_in_sub] = intersect(stim_info, info_specific{subs==s}(:,cat_stim_ids_idx),'rows','stable');
    segs(ids_in_global, (elec_ind+1):(elec_ind+elec_per_subj(subs==s)), :, :) = segs_specific{subs==s}(ids_in_sub, :, :, :);    
    elec_ind = elec_ind + elec_per_subj(subs==s);
    n_trials_per_sub(ids_in_global, subs==s) = info_specific{subs==s}.n_trials(ids_in_sub);
end
stim_info.categ = grp2idx(categorical(stim_info.categ));
stim_info.n_trials_per_sub = num2cell(n_trials_per_sub,2);

% create categories_vec - for multiple durations this is giving a unique
% number per duration x categ, otherwise it's the same as shared_stim_table.categ
switch d_str
    case 'long'
        categories_vec = stim_info.categ;
    case 'all'
        tmp = 100*unique(stim_info.categ)' + unique(stim_info.dur)/100; tmp = tmp(:)';
        cat_names = cellstr(string(cat_names) + ' ' + string(unique(stim_info.dur))); cat_names = cat_names(:); 
        categories_vec = 100*stim_info.categ + stim_info.dur/100; 
        cat_names = cat_names(ismember(tmp, unique(categories_vec)))'; categories_vec = grp2idx(categorical(categories_vec));
end

% transform n_rep == 2 to the format it's used later int the exemplar information calculations
if n_rep == 2
segs = cat(1,segs(:,:,:,1),segs(:,:,:,2));
stim_info1 = [stim_info(:,{'categ','stim_id'}), array2table(ones(size(categories_vec)),'VariableNames',{'rep'})];
stim_info2 = stim_info1; stim_info2.rep = 2*ones(size(stim_info2.rep));
stim_info = [stim_info1;stim_info2]; categories_vec = [categories_vec;categories_vec];
end

%% save

save([DATA_FOLDER, sprintf('segs_%s.mat',segs_desc)],'segs', 'stim_info', 'time_vec', 'categories_vec', 'cat_names', 'elec_info');