function [fig,ha_mat,ha_diag] = mat_diag_arrays(nrow, ncol, opt_set)
% very specific to my screen...
%
% calls Tight_subplot.m written by Pekka Kumpulainen
% (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
%
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% Written by Gal Vishne, bug reports \ requests: gal.vishne@gmail.com

if ~exist('opt_set','var'); opt_set = []; end
if ~isfield(opt_set,'bottom'); bottom = 0.35; else; bottom = opt_set.bottom; end % in mat_h units
if ~isfield(opt_set,'right'); right = 0.35; else; right = opt_set.right; end % approximate! because column widths change to match rows
if ~isfield(opt_set,'ratios'); ratios = [0.3 0.03 0.23]; else; ratios = opt_set.ratios; end % diag_h | mat_diag_gap | bet_row_gap - relative to mat_h
if ~isfield(opt_set,'w_const'); w_const = 0.02; else; w_const = opt_set.w_const; end

fig_pos = [0.05 0.05 0.09*(ncol+right) 0.225*(nrow+bottom/(1+sum(ratios)))];
fig = figure('Units','Normalized','Position',fig_pos);
set(fig,'Units','centimeters'); fig_sz = fig.Position;

mat_h = 1/(nrow*(1 + sum(ratios)) + bottom);
mat_w = (fig_sz(4)/fig_sz(3))*mat_h;
gap = [sum(ratios)*mat_h w_const]; 
marg_h =  [bottom*mat_h gap(1)];
marg_w = [(1-ncol*(mat_w+w_const)) w_const];
ha_mat = tight_subplot(nrow, ncol, gap, marg_h, marg_w);
gap(1) = (sum(ratios(2:3))+1)*mat_h;
marg_h =  [ratios(2)+bottom+1 ratios(3)]*mat_h ;
ha_diag = tight_subplot(nrow, ncol, gap, marg_h, marg_w);    
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]
end