function [CT_adj, clim] = adjust_CT(CT, clim, chance, version)
% used by plot_results_mat.m
% the idea is to get a dual sided colormap with n_colors and keep the same
% number of colors but use the range more efficiently based on clim.
% in retrospect, not clear to me why I spent time on this, but since I
% already did, enjoy...
%
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% Written by Gal Vishne, bug reports \ requests: gal.vishne@gmail.com

if ~exist('version','var')
    version = "full";
end
if round(chance-clim(1),12) == round(clim(2)-chance,12) % clims are symmetric
    CT_adj = CT;
else
clim = [min(clim(1), chance), max(clim(2), chance)];

clevels = size(CT, 1);
if version == "full"
    chance_bin = binfinder(clim, clevels, chance);
    if chance_bin == 1
        CT_adj = interpolator(CT(floor(clevels/2 + 1):clevels,:),clevels);
    elseif chance_bin == clevels
        CT_adj = interpolator(CT(1:ceil(clevels/2),:),clevels);
    else
        new_clevels = (chance_bin-1)*(clevels-chance_bin)*2; half_new_clevels = new_clevels/2;
        CT_new = interpolator(CT, new_clevels); 
        CT_adj = nan(size(CT));

        CT_adj(1:(chance_bin-1),:) = CT_new(1:half_new_clevels/(chance_bin-1):half_new_clevels,:);
        CT_adj(chance_bin,:) = CT_new(half_new_clevels,:);
        CT_adj((chance_bin+1):end,:) = CT_new(fliplr(new_clevels:(-half_new_clevels/(clevels-chance_bin)):(half_new_clevels+1)),:);
    end
elseif version == "cropped"
    clim_full = max(abs(clim-chance),[],'all');clim_full = chance+[-clim_full, clim_full];
    crop_min_or_max = find(round(clim_full,12)~=round(clim,12));
    if isempty(crop_min_or_max)
        CT_adj = CT;
    else
    % ch = chance, cl1\cl2 = clim1\clim2, nchb = new_chance_bin, nb = n_bins = clevels.
    % logic: we want to adjust cl1\cl2 (depending on the cropped side) to
    % be a bit smaller\larger so that (ch-cl1)/(cl2-ch) is exactly equal
    % to nchb/(nb-nchb). (it's not since nchb, nb is an integer and
    % ch\cl1\cl2 doesn't need to be).
    % 
    %   EQ1: (ch-cl1)/(cl2-ch) = nchb/(nb-nchb)
    % 
    % - option 1: cl1 is the cropping point -> (ch-cl1)/(cl2-ch) < 1/2, 
    % decreasing cl1 will make the nominator larger and get it closer to 1/2.
    % inverting the logic, if we take the ratio of nchb/(nb-nchb) and
    % increase it to be closer to 1/2 and then solve the equation we get the needed cl1.
    % as long as nchb is smaller than nb/2 (which is our case) it also 
    % increases when nchb is increased, so that's what we do.
    % - option 2: cl2 is the cropping point -> (ch-cl1)/(cl2-ch) > 1/2,
    % increase cl2 -> larger denominator -> smaller ratio and closer to
    % 1/2. so here we want to take the ratio of nchb/(nb-nchb) smaller by
    % decreasing nchb.
    
    % if clevels is odd the ratio should be (nchb-0.5)/(nb-nchb+0.5), which
    % changes some equations below.
    
    % how do we do this?
    % solve EQ1 for nchb and we get:
    chance_bin_unrounded = clevels * (chance-clim(1))/ range(clim);
    if isodd(clevels)
        chance_bin_unrounded = chance_bin_unrounded + 0.5;
    end
    if crop_min_or_max == 1 % crop from below
        new_chance_bin = ceil(chance_bin_unrounded); % increase nchb
        if isodd(clevels)
            new_chance_bin = new_chance_bin - 0.5;
        end
        clim(1) = (chance*clevels - new_chance_bin*clim(2))/(clevels - new_chance_bin); % and solve back for clim(1)
    else
        new_chance_bin = floor(chance_bin_unrounded); % decrease nchb
        if isodd(clevels)
            new_chance_bin = new_chance_bin - 0.5;
        end
        clim(2) = clim(1) + clevels*(chance - clim(1))/new_chance_bin; % and solve back for clim(2)
    end
    % now check if we still need to crop anything
    if round(chance-clim(1),12) == round(clim(2)-chance,12) % clims are symmetric
        CT_adj = CT;
    else
        multiplier = 1000; % to make the cut more accurate, didn't check this part in detail
        new_clevels = multiplier*clevels;
        CT_new = interpolator(CT, new_clevels); 
        crop_bin = binfinder(clim_full, new_clevels, clim(crop_min_or_max));
        if crop_min_or_max == 1 % crop from below
            ids_to_use = crop_bin:multiplier:new_clevels;
        else
            ids_to_use = fliplr((crop_bin-1):(-multiplier):1);
        end
        CT_adj = interpolator(CT_new(ids_to_use,:),clevels);

%         crop_bin = binfinder(clim_full, clevels, clim(crop_min_or_max));
%         if crop_min_or_max == 1 % crop from below
%             ids_to_use = crop_bin:clevels;
%         else
%             ids_to_use = 1:crop_bin;
%         end
%         CT_adj = interpolator(CT(ids_to_use,:),clevels);
    end
    end
end
end
CT_adj = min(max(CT_adj,0),1);
end

function bin_loc = binfinder(edges, n_bins, point)
bin_sep = linspace(edges(1), edges(2), n_bins+1);
bin_loc = find(bin_sep(2:end) > point & bin_sep(1:end-1) <= point);
end

function new_cmap = interpolator(cmap,num)
interp_num = num*max(ceil(size(cmap,1)/num),1);
new_cmap = interpolate_cbrewer(255*cmap, 'spline', interp_num)/255;
end