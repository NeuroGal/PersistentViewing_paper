function fig = mult_dur_mat_plot(res_mats, diff_mats, time_vec, res_masks, diff_masks, cluster_pvals, suptitle)
% tailored to produce figures 4 & supplementary figure 7
% Gal Vishne ~2022-23
n_dur = size(res_mats,3); 
fig = figure('Units','Normalized','Position',[0.05 0.05 0.5 0.585]);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]
gap      = [0.01 0.01];  marg_h   = [0.1 0.1];  marg_w   = [0.1 0.1];  
marg_diffs = (1-sum(marg_w)+gap(2))/(2*n_dur); marg_diffs = marg_w + [marg_diffs marg_diffs];
marg_h_1 = [(1+gap(1)-diff(marg_h))/2 marg_h(2)]; marg_h_2 = [marg_h(1) (1+gap(1)+diff(marg_h))/2];
mat_ha = tight_subplot(1, n_dur,   gap(2), marg_h_1, marg_w); diff_ha = tight_subplot(1, n_dur-1, gap(2), marg_h_2, marg_diffs); all_ax = [mat_ha;diff_ha];

clim_mat = [0 100]; clim_diff = [-50 50]; 
plot_inputs = {'decoding','no_fig', 'stim_len',[],'mask',[],'cluster_pval',[],'largest_clust'};
res_plot_inputs = [plot_inputs,{'clims',clim_mat,'chance_level',50}];
diff_plot_inputs = [plot_inputs,{'clims',clim_diff,'chance_level',0}];
for dur_i = 1:n_dur
    cur_dur = 600*(dur_i-0.5); next_dur = 600*(dur_i+0.5); keep_idx = time_vec <= next_dur;
    ticks_x = [0, cur_dur, next_dur]; 
    if dur_i == 1; ticks_y = ticks_x; else; ticks_y = next_dur; end
    adj_ratio = range(time_vec(keep_idx))/range(time_vec);
    mat_ha(dur_i).Position([3 4])  = adj_ratio*mat_ha(dur_i).Position([3 4]);
    axes(mat_ha(dur_i));
    plot_results_mat(res_mats(keep_idx,keep_idx,dur_i), time_vec(keep_idx), sprintf('%dms', cur_dur), res_plot_inputs{:},...
		'mask',res_masks(keep_idx,keep_idx,dur_i),'cluster_pval',cluster_pvals{1}{dur_i},'stim_len',cur_dur); colorbar off
    set(gca,'XTick',ticks_x,'YTick',ticks_y);
    if dur_i < n_dur
        diff_ha(dur_i).Position([3 4]) = adj_ratio*diff_ha(dur_i).Position([3 4]);
        axes(diff_ha(dur_i));
        plot_results_mat(diff_mats(keep_idx,keep_idx,dur_i), time_vec(keep_idx), sprintf('%d-%dms',next_dur, cur_dur), diff_plot_inputs{:}, ...
			'mask',diff_masks(keep_idx,keep_idx,dur_i),'cluster_pval',cluster_pvals{2}{dur_i},'stim_len',cur_dur); colorbar off
        set(gca,'XTick',ticks_x,'YTick',ticks_y);
    end
end
set(all_ax,'XLabel',[],'YLabel',[],'FontSize',12);
add_cb(mat_ha(2), 'AUC (%)'); add_cb(diff_ha(2), 'AUC diff (%)');

% final adjustments
desired_dist = 0.035;
get_w = @(ax) ax.Position(1)+ax.Position(3)+desired_dist;
mat_ha(2).Position(1) = get_w(mat_ha(1)); diff_ha(2).Position(1) = get_w(diff_ha(1));
marg_w(2) = marg_w(2) + mat_ha(3).Position(1)-get_w(mat_ha(2));
mat_ha(3).Position(1) = get_w(mat_ha(2));
mat_ha(3).Position(3) = mat_ha(2).Position(3)*mat_ha(3).Position(4)/mat_ha(2).Position(4); % no idea why this is needed, but it's messed up without this

mat_ha(2).Title.Position(2) = mat_ha(3).Title.Position(2);
mat_ha(1).Title.Position(2) = mat_ha(3).Title.Position(2);
diff_ha(1).Title.Position(2) = diff_ha(2).Title.Position(2);

supax_pos = [marg_w(1) marg_h(1) 1-sum(marg_w) 1-sum(marg_h)+ 0.04];
axes('Units','normalized','Position',supax_pos,'visible','off');
title(suptitle.txt, 'Visible','on', 'Color', suptitle.col, 'FontSize',16)
xlabel('Test time (ms)',  'Visible','on','FontSize',13);
ylabel('Train time (ms)', 'Visible','on','FontSize',13); 

function add_cb(rel_ax, tit)
cb = colorbar(rel_ax); pos = get(rel_ax,'Position');
pos([1 3]) = [0.78 0.02]; 
set(cb,'position',pos,'box','off'); title(cb,tit);
end
end