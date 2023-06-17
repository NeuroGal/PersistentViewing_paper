function fig = bar_plot_prop(prop_vals, stat_dat, txt_dat, viz_dat, prop_vals2)
% txt_dat - tit, ylab, xticks, add_legend, leg, ylab2
% stat_dat - thresh, pvals, add_stars, add_lines
% viz_dat - cmap, bar_prop, fig_w, ylims, gap, b_size_main
%   cmap after taking ROIs (n_regs x 3 x n_props) or (1 x 3 x n_props) [e.g. for categories]
%
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% Written by Gal Vishne, bug reports \ requests: gal.vishne@gmail.com

if ~isfield(viz_dat,'fig_w'); viz_dat.fig_w = 0.025; end
if ~isfield(stat_dat,'add_lines'); stat_dat.add_lines = false; end
if ~isfield(stat_dat,'add_stars'); stat_dat.add_stars = false; end
if ~isfield(txt_dat,'add_legend'); txt_dat.add_legend = false; end
if ~isfield(viz_dat,'bar_prop'); viz_dat.bar_prop = {'FaceColor','flat','LineStyle','none','FaceAlpha',0.75}; end
if exist('prop_vals2','var'); two_y_ax = true; else; two_y_ax = false;  end

star_settings = {'FontWeight','bold','HorizontalAlignment','center','FontSize',12,'Color','k'};
signif_levels = [0.05 0.01 0.001 -1];
[n_regs, n_props] = size(prop_vals);
if ~isfield(viz_dat,'gap'); gap_size = 0.08/n_props; else; gap_size = viz_dat.gap; end
if ~isfield(viz_dat,'b_size_main'); b_size_main = 0.8; else;  b_size_main = viz_dat.b_size_main; end
b_size = ((b_size_main+gap_size)/n_props) - gap_size;
tick_add = ((-b_size_main+b_size)/2):(b_size+gap_size):((b_size_main-b_size)/2);
if size(viz_dat.cmap,1) ~= 1
cmap_hsv = cellfun_wrap(@rgb2hsv, num2cell(viz_dat.cmap,2), true);
graycols = squeeze(mean(cmap_hsv(:,3,:),1));
%graycols = max(min(graycols*1.5 -0.75*mean(graycols),1),0);
viz_dat.cmap = [ones(1,3,n_props).*permute(graycols,[3 2 1]);viz_dat.cmap];
end
if ~isfield(viz_dat,'ylims') || (isfield(viz_dat,'ylims') && isempty(viz_dat.ylims))
viz_dat.ylims = limer(prop_vals,1); viz_dat.ylims(1) = 0; if stat_dat.add_stars; viz_dat.ylims(2) = viz_dat.ylims(2)+1; end
end
if two_y_ax && (~isfield(viz_dat,'ylims2') || (isfield(viz_dat,'ylims2') && isempty(viz_dat.ylims2)))
viz_dat.ylims2 = limer(prop_vals2,1/10);
end
fig = figure('Units','Normalized','Position',[0.2 0.2 0.03+viz_dat.fig_w*n_props*n_regs 0.3]);
for prop = 1:n_props
cur_x = (1:n_regs) + tick_add(prop);
xplt = cur_x; yplt = prop_vals(:,prop);
if size(viz_dat.cmap,1)~=1; xplt = [-2,xplt]; yplt = [0;yplt]; end
if two_y_ax; yyaxis left; end
b_h(prop) = bar(xplt,yplt,b_size,'CData',viz_dat.cmap(:,:,prop),viz_dat.bar_prop{:});hold on
if stat_dat.add_lines
for loc_i = 1:n_regs
    plt = plot(cur_x(loc_i)+[-1,1]*(b_size+max(gap_size-0.01,0))/2,[stat_dat.thresh(loc_i,prop) stat_dat.thresh(loc_i,prop)], 'Color', [1 1 1]*0.6, 'LineWidth', 1.5);
end
end
if stat_dat.add_stars
for sig = 1:3
    signif_idx = stat_dat.pvals(:,prop)<=signif_levels(sig) &  stat_dat.pvals(:,prop)>signif_levels(sig+1);
    text(cur_x(signif_idx),(viz_dat.ylims(2)-0.05*range(viz_dat.ylims))*ones(1,sum(signif_idx)),repmat('*',1,sig),star_settings{:});
end
end
if two_y_ax; yyaxis right;
cmap = viz_dat.cmap(:,:,prop); if size(cmap,1) ~= 1; cmap = cmap(2:end,:); end
scatter(cur_x,prop_vals2(:,prop),70,cmap,'filled')
end
end
if two_y_ax
yyaxis right
ylabel(txt_dat.ylab2); ylim(viz_dat.ylims2);
ha = gca; set(ha.YAxis,'Color','k')    
yyaxis left
end
title(txt_dat.tit); ylabel(txt_dat.ylab); xlim([0.4 n_regs+0.6]); ylim(viz_dat.ylims);
set(gca,'XTick',1:n_regs,'XTickLabels',txt_dat.xticks,'FontSize',12.5,'Box','off')
if txt_dat.add_legend
    if n_props == 1; leg_plt = plt; else; leg_plt = b_h; end
    legend(leg_plt,txt_dat.leg,'FontSize',10,'Box','off');
end
end