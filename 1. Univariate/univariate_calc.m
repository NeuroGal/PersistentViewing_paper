%% Analyses for results section 1 (Figure 1 & Supp Figure 1)
% take the output from here to univariate_fig.m
% comparisons between regions and summary stats are also found there

clear all;clc;close all

%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));
measures = {'peak-time','peak-amp','attenuation'};
types    = {'Responsive','Selective'};

%% Responsiveness segments
% create mean segments across categories (weighted according to the number of exemplars in each categories,
% only including categories each electrode was responsive for)

% Load segments created with specific_segs (d_str corresponds to tv_end = 900 there)
n_rep           = 1;
d_str           = 'long';
subs            = 1:10;
sample_hz       = 1000;
cat_nums        = 1:4;
av_in_cat       = true;
segs_desc = sprintf('%drep_dur%s_subs%s_%dHz_cat%s_av%s',n_rep,d_str,join(string(subs),''),...
    sample_hz,join(string(cat_nums),''),string(av_in_cat));
load([DATA_FOLDER, sprintf('segs_%s.mat',segs_desc)]);
n_resp_e = sum(elec_info.IsResponsive);

tmp = cell2mat(stim_info.n_trials_per_sub);
mult_nums = tmp(:,elec_info.Patient);
mult_nums = mult_nums.*elec_info.RespPerCat_FWOA'; % just average what it is responsive for
segs_resp = squeeze(sum(mult_nums.*segs,1)./sum(mult_nums,1)); % elec x time
segs_resp = segs_resp(elec_info.IsResponsive,:);

%% Selectivity segments
% for this we load the full basic segments because we want to use all trials
% possible for each patients not just the ones shared across patients (in
% the responsiveness segments we do the averaging in 'specific_segs')
% ! do load the specific segs before this (for time_vec etc)

if ~exist('segs_basic','var'); load([DATA_FOLDER, sprintf('basic_segs_%dHz.mat',sample_hz)],'segs_basic','trial_info_basic','time_vec_basic'); end
segs_select  = nan(n_resp_e, length(time_vec));

subj_elec_idx = 0;
for s = subs
    time_filt = time_vec_basic > -100 & time_vec_basic <= 900; % change the output time range
    trial_info = trial_info_basic{s}; dur_filt = trial_info.durations >= 900;  % select relevant durations
    e_filt = elec_info.IsResponsive(elec_info.Patient==s);
    segs = double(segs_basic{s}(dur_filt,e_filt,time_filt)); 
    cat_nums = categorical(trial_info.cat_nums(dur_filt));
    n_elec = size(segs, 2);
    for t = 1:length(time_vec)
        for e = 1:n_elec
            [~,tmp] = anova1(segs(:, e, t),cat_nums,'off');
            segs_select(subj_elec_idx+e, t) = 100*tmp{2,2}/tmp{4,2};
        end
    end
    subj_elec_idx = subj_elec_idx + n_elec;
    fprintf('Done with subj #%d\n',s)
end

%% calculate univariate properties (peak time, peak value, attenuation)

univar_prop = nan(n_resp_e, 3, 2); % elec, meas(peak time, peak amp, attenuation), resp\select
for type = 1:2
    if type == 1
        signs = [1,-1]; segs = segs_resp;
    elseif type == 2
        signs = 1; segs = segs_select;
    end
    for e_sign = signs
        if type == 1; e_filt = elec_info.RespSign(elec_info.IsResponsive) == e_sign; else; e_filt = true(sum(elec_info.IsResponsive),1); end
        segs_tmp = segs(e_filt, time_vec>0);
        [peak_val, peak_idx] = max(e_sign.*segs_tmp, [], 2);
        peak_val = e_sign.*peak_val;
        univar_prop(e_filt, 1, type) = time_vec(find(time_vec == 0) + peak_idx);
        univar_prop(e_filt, 2, type) = peak_val;
        segs_end = mean(segs(e_filt, time_vec>=800 & time_vec<=900),2);
        univar_prop(e_filt, 3, type) = 100*(peak_val-segs_end)./peak_val;
    end
end

%% save (for use in univariate_fig)

elec_info = elec_info(elec_info.IsResponsive,:);
save([DATA_FOLDER,'univariate_calc.mat'],'segs_resp','segs_select','univar_prop','measures','types','elec_info','time_vec');
