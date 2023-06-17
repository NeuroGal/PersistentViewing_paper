function [corr_r, mask, corr_p, cluster_p, corr_perm_r] = corr_w_cat_mod(rdm, model_rdm, stat_info, corr_type)
% Compute correlation with category models as done in this paper:
%   Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% Please cite if used.
%
% Inputs:
%   Mandatory:
%     rdm       - n_stim x n_stim x n_time
%     model_rdm - n_stim x n_stim
%     stat_info - see rdm_corr_stats.m for more details, here we need:
%               stat_type - used to decide if we extract permutation data
%               n_perm - how many permutations to run
%               n_sides - 1/2-sided test (1-sided is only right-sided)
%   Optional:
%     corr_type - string indicating which correlation measure to use. default: Spearman. 
%
% Written by Gal Vishne, lab of Leon Y. Deouell, ~2022
% Send bug reports and requests to gal.vishne@gmail.com

if ~exist('corr_type','var'); corr_type = 'Spearman'; end
[stat_info, tail_type] = parse_stat_info(stat_info);
corr_inputs = {'Type',corr_type,'Tail',tail_type};
rdm = reshape_dists(rdm); % turn from n_stim x n_stim x n_time to n_dists x n_time (keeping only upper triangle)
[corr_r, corr_p] = get_corr_vals(rdm, model_rdm, corr_inputs, false); % last input is dont permute

if ~strcmp(stat_info.stat_type,'analytic')
    dims = [size(corr_r), stat_info.n_perm]; corr_perm_r = nan(dims); corr_perm_p = nan(dims);
    for p = 1:stat_info.n_perm
        [corr_perm_r(:,:,p), corr_perm_p(:,:,p)] = get_corr_vals(rdm, model_rdm, corr_inputs, true); % last input is yes to permute
        if mod(p,100)==0; fprintf('Done with permutation %d/%d\n',p,stat_info.n_perm); end
    end
else
    corr_perm_r = []; corr_perm_p = []; 
end
[mask, corr_p, cluster_p] = rdm_corr_stats(stat_info, corr_r, corr_perm_r, corr_p, corr_perm_p);
end

function [corr_r, corr_p] = get_corr_vals(rdm, model_rdm, corr_inputs, perm)
perm_n = size(model_rdm,2); 
if perm; perm_order = randperm(perm_n); else; perm_order = 1:perm_n; end
model_rdm = reshape_dists(model_rdm(perm_order,perm_order,:)); % rdm reshaped outside
[corr_r, corr_p] = corr(rdm, model_rdm, corr_inputs{:});
end

function [stat_info, tail_type] = parse_stat_info(stat_info)
stat_info.get_diag_stats = false; stat_info.chance_level = 0; % adding relevant fields
fields = {'stat_type','n_sides','n_perm'};
defaults = {'analytic', 1, 1000}; 
for f = 1:length(fields)
    if ~isfield(stat_info, fields{f}) 
        stat_info.(fields{f}) = defaults{f};
        warning('Setting %s to %s', fields{f}, string(defaults{f}))
    end
end
if stat_info.n_sides == 1
    tail_type = 'right';
elseif stat_info.n_sides == 2
    tail_type = 'both';
end
end