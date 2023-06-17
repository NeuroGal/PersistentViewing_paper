function fig = dec_diag_plot(reg_dat, time_vec, reg_info, comp_info, subj_mode, fig_inputs)
% tailored for plotting things I needed for Supp 6
% reg_info - names, ord (i.e. indices to plot and which order), reg_cmap
% comp_info - names, ord (same)
% subj_mode - is this single subject data plotted for one region across subjects?

if ~exist('subj_mode','var'); subj_mode = false; end
all_diags = cell2mat(reg_dat.results(comp_info.ord,reg_info.ord,2));
ylims = limer(all_diags, 1/5); ylims(2) = 100;
ncol = length(comp_info.ord); nrow = length(reg_info.ord); fig_pos = [0.05 0.05 0.09*(2+ncol) 0.08*(1+nrow)];
gap = [0.008 0.005]; marg_h = [0.15 0.1]; marg_w = [0.5 0.05]; marg_w = marg_w/(sum(marg_w)*(1+ncol)); marg_h = marg_h/(sum(marg_h)*(1+nrow));
[fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w);
tick_ax_num = ncol*(nrow-1)+1;
for reg_i = reg_info.ord
    if ~subj_mode; col = reg_info.reg_cmap(reg_i,:); else; col = reg_info.reg_cmap; end
    for comp_i = comp_info.ord
        plt_idx = (find(reg_info.ord==reg_i)-1)*ncol+find(comp_info.ord==comp_i);
        vals = reg_dat.results{comp_i,reg_i,2}; mask = reg_dat.masks_clust{comp_i,reg_i,2}; pval = reg_dat.clust_ps{comp_i,reg_i,2}; % assumes the calculation code already ran cluster based statistics
        axes(ha(plt_idx));
        nice_line_plot(vals, time_vec, '', 'color',col, fig_inputs{:}, 'ylims',ylims,'cluster_pval',pval,'mask',mask);
        if reg_i == reg_info.ord(1);title(comp_info.names{comp_i},'Color','k'); end
        if plt_idx~=tick_ax_num; set(gca,'XLabel',[],'YLabel',[],'XTick',[],'YTick',[]); end
        if comp_i == comp_info.ord(2)
            if ~subj_mode
            ylb = ylabel(split(reg_info.names(reg_i)),'Color',col);
            else
            ylb = ylabel("S"+string(reg_i));
            end
            set(ylb,'Rotation',0,'FontWeight','bold');
            ylb.Position(1) = -1500-25*ncol;
        end
    end
end
end