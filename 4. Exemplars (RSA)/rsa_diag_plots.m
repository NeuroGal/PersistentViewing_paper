function figs = rsa_diag_plots(ROIs, RELs, diag_info, viz_dat, fig_title_func)
% diag_info: models x reliability x region, has all of the fields of reg_dat.stats
% (inserting category data to this function is done with one category per region so it's not in the models axis)
% viz_dat: plot_inputs (gen_plot_inputs_line), ROInames_full, font_size, time_vec
%   ylims (optional; first dim for reliability metric, i.e. ylims(rel_i,:))
%   cmap (reg_cmap\ cmap with models in 3rd arg), tit_cmap (for the case with models, probably reg_cmap)
nrow = 1; if length(ROIs) > 5; nrow = 2; end; ncol = ceil(length(ROIs)/nrow);
fig_pos = [0.1 0.25 0.025+0.155*ncol 0.08+0.13*nrow];
marg_h = [0.05+0.22/nrow 0.015+0.12/nrow]; marg_w = [0.05+0.06/ncol 0.01+0.03/ncol]; gap = [0.105+0.08/nrow 0.02];
figs = []; n_model = size(diag_info,1);
for rel_i = RELs
    [fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w);
    all_diags = cellfun_wrap(@(x) x.vals, permute(diag_info(:,rel_i,ROIs),[2 1 3]),true); % time x models x region
    ylims = limer(all_diags,1);
    if isfield(viz_dat,'ylims'); ylims = viz_dat.ylims(rel_i,:); end
    plot_inputs = [viz_dat.plot_inputs,{'ylims',ylims}];
    for reg_i = ROIs
        axes(ha(ROIs==reg_i));
        for mod_i = 1:n_model
            cur_idx = sub2ind(size(diag_info),mod_i,rel_i,reg_i);
            rel_met = diag_info{cur_idx}.vals; mask = diag_info{cur_idx}.mask; cluster_p = diag_info{cur_idx}.cluster_p;
            if n_model==1; added_inputs = {'mask', mask, 'cluster_pval', cluster_p}; else; added_inputs = {}; end
            nice_line_plot(rel_met, viz_dat.time_vec, viz_dat.ROInames_full(reg_i), plot_inputs{:}, added_inputs{:}, 'color', viz_dat.cmap(reg_i,:,mod_i)); hold on
            if n_model>1; add_mask_stars(mask, cluster_p, rel_met, viz_dat.time_vec, ylims, viz_dat.cmap(reg_i,:,mod_i), [50/range(viz_dat.time_vec) 0.035*mod_i], 'left', 11); end
        end
        if reg_i ~= ROIs(1); set(gca,'Xlabel',[],'Ylabel',[],'YTick',[]); end
        if n_model>1; title(viz_dat.ROInames_full(reg_i),'Color',viz_dat.tit_cmap(reg_i,:)); end
    end
    set(ha,'FontSize',viz_dat.font_size);set(fig,'Name', fig_title_func(rel_i)); figs = [figs,fig];
end
end