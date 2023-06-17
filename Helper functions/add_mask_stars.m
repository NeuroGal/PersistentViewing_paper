function add_mask_stars(mask, pval, vals, tv, ylims, col, offset, hor_align, fontsize)
% add lines for each cluster & adjacent stars marking the p value (corresponding to the temporal order)
%
% pval - cluster p
% vals -!!! of the stat used to define the cluster! so it will arrange the stars based on cluster size
% offset (x,y - from the left & bottom) in normalized units
%
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% Written by Gal Vishne, bug reports \ requests: gal.vishne@gmail.com

if ~exist('hor_align', 'var'); hor_align = 'center'; end
if ~exist('fontsize', 'var'); fontsize = 13; end

mask_plt = nan(size(mask)); mask(isnan(mask)) = false; mask_plt(boolean(mask)) = true;
plot(tv, mask_plt*(ylims(1)+offset(2)*range(ylims)), 'LineWidth', 2, 'Color', col);
signif_levels = [0.05 0.01 0.001 0];
if ~isempty(pval)
    if numel(pval)>1
        pval = sort(pval); % so it will be from largest cluster to smallest cluster, and then the ordering makes senses
        tmp = bwconncomp(mask); [~,sort_idx] = sort(cellfun(@(x) sum(vals(x)), tmp.PixelIdxList),'descen'); % resort pstr according to the order of the clusters
        pval = pval(sort_idx); pval(pval>=0.05)=[];
        pstr = join(string(arrayfun(@(x) repmat('*',1,find(x<=signif_levels,1,'last')),pval,'UniformOutput',false)),' ');
    else
        pstr = repmat('*',1,find(pval<=signif_levels,1,'last'));
    end
else
    pstr = '';
end
u_change = @(offset,lims) lims(1)+offset*range(lims);
text(u_change(offset(1),tv([1 end])),u_change(offset(2)-0.005,ylims),pstr,'FontSize',fontsize,'FontWeight','bold','HorizontalAlignment',hor_align,'Color',col);
end