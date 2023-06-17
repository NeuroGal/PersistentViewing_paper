function figs = rsa_mat_plots(ROIs, RELs, reg_dat, viz_dat, fig_title_func, use_clust)
% viz_dat: plot_inputs (gen_plot_inputs_mat), time_vec, ROInames_full, metric, cmap (reg_cmap\ cmap with models)
%   clims (optional; first dim for reliability metric, i.e. clims(rel_i,:))
% for figure 10d-e this function is used per ROI across models, so ROIs is actually model indices, the dimensions of reg_dat are flipped accordingly etc. 
if ~exist('use_clust','var'); use_clust = false; end
nrow = 1; if length(ROIs) > 5; nrow = 2; end; ncol = ceil(length(ROIs)/nrow);
fig_pos = [0.1 0.25 0.025+0.155*ncol 1.5*(0.08+0.13*nrow)];
marg_h = 0.75*[0.05+0.22/nrow 0.015+0.12/nrow]; marg_w = [0.1+0.12/ncol 0.01+0.03/ncol]; gap = [0.105+0.08/nrow 0.02];
figs = []; 
for rel_i = RELs
    [fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w);
    all_mats = reg_dat.relis(:,:,:,rel_i,ROIs);
    clims = limer(all_mats,1);
    if isfield(viz_dat,'clims'); clims = viz_dat.clims(rel_i,:); end
    plot_inputs = [viz_dat.plot_inputs,{'no_fig','clims',clims}];
    for reg_i = ROIs
        axes(ha(ROIs==reg_i));
        rel_met = all_mats(:,:,:,:,ROIs==reg_i);
        if ~use_clust
            added_inputs = {'mask', reg_dat.masks_perm(:,:,:,rel_i,reg_i)};
        else
            added_inputs = {'mask', reg_dat.masks(:,:,:,rel_i,reg_i), 'cluster_pval', reg_dat.clust_ps{:,rel_i,reg_i}};
        end
        plot_results_mat(rel_met, viz_dat.time_vec, viz_dat.ROInames_full(reg_i), added_inputs{:}, plot_inputs{:}, 'line_col', viz_dat.cmap(reg_i,:)); colorbar off
        if find(ROIs==reg_i)~=ncol*(nrow-1)+1
            set(gca,'Xlabel',[],'Ylabel',[],'YTick',[],'YLabel',[]);
        else
            if rel_i == 1
                ylabel('Time same repetition (ms)'); xlabel('Time other repetiton (ms)'); 
            elseif rel_i == 2
                ylabel('Time Geometry 1 (ms)'); xlabel('Time Geometry 2 (ms)'); 
            end
        end
    end
    cb = colorbar(ha(1)); ax_pos = ha(1).Position; set(cb, 'Position', [0.0425 ax_pos(2)+ax_pos(4)*0.1 0.025-ncol*0.0025 ax_pos(4)*0.8],'Box','off'); title(cb,viz_dat.metric);
    set(fig,'Name', fig_title_func(rel_i)); figs = [figs,fig];
end
end