function [fig,ha] = fig_ax_wrap(fig_pos, nrow, ncol, gap, marg_h, marg_w)
% calls Tight_subplot.m written by Pekka Kumpulainen
% (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
fig = figure('Units','Normalized','Position',fig_pos);
ha = tight_subplot(nrow, ncol, gap, marg_h, marg_w);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]
end
