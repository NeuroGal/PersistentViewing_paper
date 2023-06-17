function fig = dec_mat_diag_plot(reg_dat, time_vec, reg_info, comp_info, fig_info, fig_inputs)
% tailored for plotting things I needed for Supp 5-6
% reg_info - names, ord (i.e. indices to plot and which order), reg_cmap
% comp_info - names, ord (same)
% fig_info - nrow, ncol, opt_set (see mat_diag_arrays), del_idx (specific for how I wanted to arrange plots)

comp_mode = length(comp_info.ord)>1; % true -> we do this as reg x comp, otherwise it's just each reg (or subj) one comp
if ~isfield(fig_info,'del_idx'); fig_info.del_idx = []; end

[fig,ha_mat,ha_diag] = mat_diag_arrays(fig_info.nrow, fig_info.ncol, fig_info.opt_set);
delete([ha_mat(fig_info.del_idx);ha_diag(fig_info.del_idx)]); ha_mat(fig_info.del_idx) = []; ha_diag(fig_info.del_idx) = []; 
tick_ax_num = fig_info.ncol*(fig_info.nrow-1)+1; tick_ax_num = tick_ax_num - sum(fig_info.del_idx<tick_ax_num);

all_diags = [reg_dat.results{comp_info.ord,reg_info.ord,2}];
ylims = limer(all_diags, 1/5); clims = [0 100];
 
for reg_i = reg_info.ord
    for comp_i = comp_info.ord
        if comp_mode
            plt_idx = (find(reg_info.ord==reg_i)-1)*length(comp_info.ord)+find(comp_info.ord==comp_i); col = reg_info.reg_cmap(reg_i,:);
        else
            plt_idx = find(reg_info.ord==reg_i); col = reg_info.reg_cmap;
        end
        
        axes(ha_mat(plt_idx));
        plot_results_mat(reg_dat.results{comp_i,reg_i,1}, time_vec, '', fig_inputs.mat{:},'mask',reg_dat.masks_perm{comp_i,reg_i,1},'clims',clims,'no_fig'); colorbar off
        if plt_idx~=tick_ax_num;set(gca,'XTick',[],'YTick',[],'XLabel',[],'YLabel',[]);end   
        if comp_mode && comp_i == comp_info.ord(2)
            ylb = ylabel(split(reg_info.names(reg_i)),'Color',col);
            set(ylb,'Rotation',0,'FontWeight','bold'); ylb.Position(1) = min(time_vec)-1.85*range(time_vec);
        end        

        axes(ha_diag(plt_idx));
        nice_line_plot(reg_dat.results{comp_i,reg_i,2}, time_vec,'','mask',reg_dat.masks_clust{comp_i,reg_i,2},'cluster_pval',reg_dat.clust_ps{comp_i,reg_i,2},...
            'ylims',ylims, fig_inputs.line{:},'color', col);
        set(gca,'XTick',[],'XLabel',[]);
        if plt_idx~=tick_ax_num; set(gca,'YTick',[],'YLabel',[]); else; set(gca,'YTick',[50 100]); end
        if comp_mode
            if reg_i == reg_info.ord(1); title(comp_info.names(comp_i),'Color','k'); end
        else
            title(reg_info.names(plt_idx),'Color','k');
        end
    end
end
cb = colorbar(ha_mat(1)); ax_pos = ha_mat(1).Position;
set(cb, 'Position', [0.0425 ax_pos(2)-ax_pos(4)*0.55 0.03-fig_info.ncol*0.003 ax_pos(4)*0.9],'Box','off'); title(cb,'AUC (%)');
end