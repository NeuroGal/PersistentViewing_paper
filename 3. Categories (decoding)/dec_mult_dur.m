function data_out = dec_mult_dur(data, time_vec, cat_inputs, add_inf)
% Input: data        - n_stim x n_elec x n_time (\x2)
%        time_vec    - corresponding to n_time
%        cat_inputs  - {'categories', categories_vec, cat_names}
%        add_inf     - struct with: 
%                       cfg_decoding
%                       cfg_stats
%                       mask_data
%
% This code was written as part of the analysis of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
%
% Written by Gal Vishne, lab of Leon Y. Deouell, ~2021-2022
% Send bug reports and requests to gal.vishne@gmail.com

n_time = length(time_vec);
categories_vec = cat_inputs{2}; cat_names = cat_inputs{3};
if size(cat_names,2)>size(cat_names,1)
    cat_names = cat_names';
end
tmp_cat_names = split(cat_names);
[dur_unique,~,unique_idx_dur] = unique(tmp_cat_names(:,2), 'stable');
cat_unique = unique(tmp_cat_names(:,1), 'stable');
n_dur = length(dur_unique); dur_unique = str2double(dur_unique);
n_cat = length(cat_unique);

add_stats = add_inf.mask_data; cfg_decoding = add_inf.cfg_decoding; cfg_stats = add_inf.cfg_stats;
% make the mvpa light default explicit
if ~isfield(cfg_decoding, 'metric'); cfg_decoding.metric = 'acc'; end
if ~isfield(cfg_decoding, 'classifier'); cfg_decoding.classifier = 'lda'; end
if iscell(cfg_decoding.metric)
    warning('Using only the first requested metric %s', cfg_decoding.metric{1})
    cfg_decoding.metric = cfg_decoding.metric{1};
end
if strcmp(cfg_decoding.metric,'confusion')
    cfg_decoding.metric = 'acc'; % notice the forced change
    warning('Confusion matrix is not possible. Changed metric to accuracy.')
end
if n_cat > 2 && ~strcmp(cfg_decoding.classifier,'multiclass_lda')
    cfg_decoding.classifier = 'multiclass_lda'; % notice the forced change
    warning('Changed classifier to multiclass_lda')
end
if strcmp(cfg_decoding.metric,'acc'); cfg_stats.acc = true; else; cfg_stats.acc = false; end
cfg_stats_diff = cfg_stats; cfg_stats_diff.acc = false; cfg_stats_diff.chance_level = 0; % even if the metric is acc the difference shouldn't be treated as such

n_perm = cfg_stats.n_perm;
% this is the actual action
dims = [n_time, n_time, n_dur];
res_mats = nan(dims); res_masks = nan(dims); pvals_res = cell(n_dur,2); % mask\diag
res_diag_masks = nan(n_time, 1, n_dur);
% for diff stats:
perms = nan(n_time, n_time, n_perm, n_dur);
for dur = 1:n_dur
    max_time = dur_unique(dur)+600;
    rel_stim = ismember(categories_vec,find(unique_idx_dur==dur)); keep_idx = time_vec <= max_time;
    dec_inputs = {data(rel_stim, :, keep_idx, :), categories_vec(rel_stim), cfg_decoding};
    if add_stats; dec_inputs = [dec_inputs, cfg_stats]; end
    [result_tmp, stat_tmp, ~, perm_tmp] = decoding_wrapper(dec_inputs{:}); % mask tmp is all the stats
    if iscell(result_tmp); result_tmp = result_tmp{:}; end
    res_mats(keep_idx, keep_idx, dur) = result_tmp(keep_idx,keep_idx);
    perm_tmp = perm_tmp{:}; if ~isempty(perm_tmp); perms(keep_idx, keep_idx, :, dur) = perm_tmp; end
    if add_stats
    if strcmp(cfg_decoding.classifier,'multiclass_lda')
        pvals_res(dur,:) = stat_tmp(2,:,:);
        res_masks(keep_idx, keep_idx, dur) = stat_tmp{1,:,1};
        res_diag_masks(keep_idx, 1, dur) = stat_tmp{1,:,2};
    else
        pvals_res(dur,:) = stat_tmp(2,:,:,:);
        res_masks(keep_idx, keep_idx, dur) = stat_tmp{1,:,:,1};
        res_diag_masks(keep_idx, 1, dur) = stat_tmp{1,:,:,2};
    end
    end
end
diff_masks = nan(n_time, n_time, n_dur-1); diff_diag_masks = nan(n_time, 1, n_dur-1);
pvals_diff = cell(n_dur-1,2); % mat\diag
if add_stats % calculate the diff stats
    for dur = 2:n_dur
        keep_idx = time_vec <= dur_unique(dur-1)+600; % shorter determines the extent
        vals = diff(res_mats(keep_idx,keep_idx,(dur-1):dur),[],3);
        perms_tmp = squeeze(diff(perms(keep_idx,keep_idx,:,(dur-1):dur),[],4));
        diag_vals = diag(vals); diag_perms = get_diags(perms_tmp);
        [diff_masks(keep_idx, keep_idx, dur-1), pvals_diff{dur-1,1}] = decoding_stats(cfg_stats_diff, vals, perms_tmp);
        [diff_diag_masks(keep_idx,1,dur-1), pvals_diff{dur-1,2}]   = decoding_stats(cfg_stats_diff, diag_vals, diag_perms);
    end
end
data_out.dur_unique      = dur_unique;
data_out.res_mats        = res_mats; % diff_mats & diags can easily be computed from this
data_out.res_masks       = res_masks;
data_out.diff_masks      = diff_masks;
data_out.pvals_res       = pvals_res;
data_out.pvals_diff      = pvals_diff;
data_out.res_diag_masks  = res_diag_masks;
data_out.diff_diag_masks = diff_diag_masks;
data_out.perms           = perms;
data_out.cfg_stats_diff  = cfg_stats_diff;
end