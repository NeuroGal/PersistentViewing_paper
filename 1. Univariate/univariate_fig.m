%% Figures for results section 1 (Figure 1 & Supp Figure 1)
% and also summary stats and comparisons between regions

clear all;clc;close all
%DATA_FOLDER = pwd; % change to the folder where you put basic_segs.mat and the other data
%CODE_FOLDER = pwd; % change to the folder with all of the github code (from all toolboxes)
%addpath(genpath(CODE_FOLDER));
load([DATA_FOLDER,'univariate_calc.mat'])
load([DATA_FOLDER,'figure_settings.mat'])

%% specific settings you wish to plot

type    = 1; % 1-responsiveness, 2-selectivity
meas    = 1; % 1-peak time, 2-peak amp, 3-attenuation - matters only for the violin plots
e_sign  = 1; % 1\-1, only relevant for responsiveness
ROIs    = 1:4; % in this order:  Occ, VT, Par, PFC, SM, LT, Occ-retin, Occ-not-retin, LPFC, OFC
% note that if you compare regions 7-10 you can't also show Occ & PFC
if type == 1
    rel_segs = segs_resp; e_filt = elec_info.RespSign == e_sign;
elseif type == 2
    rel_segs = segs_select; e_filt = elec_info.IsSelective;
end
if all(ROIs<=6); ROI_idx = grp2idx(elec_info.ROI); else; ROI_idx = grp2idx(elec_info.ROIsplit); end % use ROIsplit if any of regions 7-10 requested
e_filt = e_filt & ismember(ROI_idx,ROIs);
rel_prop = univar_prop(e_filt, meas, type); rel_segs = rel_segs(e_filt,:);
ROI_idx = grp2idx(categorical(ROI_idx(e_filt))); reg_names = ROInames(ROIs);

fig_settings = [];
fig_settings.reg_cmap = reg_cmap(ROIs,:);
fig_settings.ROInames = ROInames_full(ROIs);
if type == 1; type_txt = 'Response'; elseif type == 2; type_txt = 'Selectivity'; end

%% just stats, no plot yet

clc
% mean,sem,count per region (& for all regions together)
tbl_per_reg = table(ROI_idx, rel_prop, 'VariableNames',{'loc','val'});
tbl_all = table(11*ones(size(ROI_idx)), rel_prop, 'VariableNames',{'loc','val'});
summary_stats = [grpstats(tbl_per_reg,'loc',{'mean','sem'});grpstats(tbl_all,'loc',{'mean','sem'})];
summary_stats.Row = [reg_names,"ALL"];
summary_stats(:,2:end)

% anova + pairwise comparisons
[p,tbl,stats] = anova1(rel_prop,ROI_idx,'off'); % \kruskalwallis also possible
fprintf('ANOVA results: F(%d,%d) = %0.2f, p < %0.4g\n', tbl{2,3}, tbl{3,3}, tbl{2,5}, tbl{2,6})
mlt_cmp = multcompare(stats,'Display','off');
p_array = mlt_cmp(:,6);
cohens_d_array = cellfun(@(rois) cohens_d(rel_prop(ROI_idx == rois(1)),rel_prop(ROI_idx == rois(2))),num2cell(mlt_cmp(:,1:2),2));
str_ver = string(arrayfun(@(x) sprintf('%0.4g',x),[p_array,cohens_d_array],'UniformOutput',false));
comp_text_array = join(reg_names(mlt_cmp(:,1:2)),'-')+": p < " +str_ver(:,1) + ", Cohen's d = "+str_ver(:,2);
fprintf('Pairwise comparisons: \n%s\n', join(comp_text_array,newline))

%% Violin plots 
% -> run anova (above) first to have the state for the figure
% Settings per panel in the paper: (set above)
%   ROIs = 1:4 in all cases
%   Figure 1d:  type = 1; meas = 3; e_sign = 1;
%   Figure 1f:  type = 2; meas = 3;
%   Figure S1a: type = 1; e_sign = 1; meas = 1; 
%   Figure S1a: type = 1; e_sign = 1; meas = 2; 
%   Figure S1g: type = 2; meas = 1; 
%   Figure S1h: type = 2; meas = 2; 

fig_name = sprintf('Univar %s %s %s',types{type}, measures{meas}, join(reg_names));
if type == 1; fig_name = sprintf('%s sign%s',fig_name,string(e_sign)); end
fig_settings.fig_name = fig_name;

violin_titles = {'Peak %s time', 'Peak %s magnitude', '%s attenuation'};
if meas ~= 3; type_txt = lower(type_txt); end
fig_settings.title = sprintf(violin_titles{meas},type_txt);

if meas == 1
    ylab = 'ms';
elseif meas == 2
    if type == 1; ylab = 'HFA (a.u.)'; else; ylab = '% explained variance'; end
elseif meas == 3
    ylab = '% from peak';
end
fig_settings.ylab = ylab;

fig = violin_panel(rel_prop, ROI_idx, fig_settings, mlt_cmp);

%% Dynamics plot
% Settings per panel in the paper:
%   Figure 1e:   type = 1; e_sign = 1; ROIs = 1:4
%   Figure 1g:   type = 2; ROIs = 1:4;
%   Figure S1c:  type = 1; e_sign = 1; ROIs = 5:6;
%   Figure S1d:  type = 1; e_sign = 1; ROIs = 7:10;
%   Figure S1e:  type = 1; e_sign = -1; ROIs = 1:6;

dynamics_ylabs = {'HFA amplitude (a.u.)','Category selectivity index (% explained variance)'};
fig_settings.title = sprintf('%s dynamics',type_txt);
fig_settings.fig_name = sprintf('Univar %s dynamics %s',types{type},join(reg_names));
fig_settings.ylab = dynamics_ylabs{type};

fig = dynamics_panel(rel_segs, time_vec, ROI_idx, fig_settings);

%% functions

function fig = violin_panel(rel_prop, ROI_idx, settings, mlt_cmp)
% settings needs to include: ylab, title, fig_name, ROInames, reg_cmap
n_reg = length(unique(ROI_idx));
settings.ylim = limer(rel_prop); % see helper_functions folder

fig = figure('Units','Normalized','Position',[0.2 0.3 (n_reg+1.1)*0.038 0.4],'Name',settings.fig_name);
violin_array = violinplot(rel_prop,ROI_idx,'ShowMean',true); hold on
for roi_i = 1:n_reg; violin_array(roi_i).ViolinColor = settings.reg_cmap(roi_i,:); end
ylim(settings.ylim); xlim([0.5 n_reg+0.5]);
ylabel(settings.ylab); title(settings.title);
set(gca,'XTick',1:n_reg,'XTickLabels',settings.ROInames,'FontSize',15)

% add stats lines & stars
signif = find(mlt_cmp(:,6)<=0.05)';
for sig = signif
    n_stars = sum(mlt_cmp(sig,6)<=[0.05 0.01 0.001]);
    yrange = range(settings.ylim);
    y_pos = settings.ylim(2)-sig*0.05*yrange;
    plot(mlt_cmp(sig,1:2)+[0.1 -0.1],[y_pos y_pos],'k')
    text(mean(mlt_cmp(sig,1:2)),y_pos+0.01*yrange,repmat('*',1,n_stars),'HorizontalAlignment','center','FontSize',14);
end
end

function fig = dynamics_panel(rel_segs, time_vec, ROI_idx, settings)
% settings needs to include: ylab, title, fig_name, ROInames, reg_cmap
n_reg = length(unique(ROI_idx));

% get ylimits for the main plot
mean_vals = cellfun_wrap(@(roi_i) mean(rel_segs(ROI_idx==roi_i, :),1),num2cell((1:n_reg)'),true);
ylimit = limer(mean_vals*1.1); % 1.1 to account for the SEM shaded area

% plot dynamics
fig = figure('Units','Normalized','Position',[0.3 0.3 0.35 0.35],'Name',settings.fig_name);
plot([0 0], ylimit, 'k:'); hold on
if contains(settings.ylab, 'HFA'); plot([time_vec(1) time_vec(end)],[0 0],'k:'); end
e_nums = nan(1,n_reg);
for roi_i = 1:n_reg
    inputs = {time_vec, rel_segs(ROI_idx==roi_i, :)', 'Color', settings.reg_cmap(roi_i,:),'LineWidth',1.5};
    if isvector(inputs{2}); plot(inputs{:}); else; varplot(inputs{:}); end
    e_nums(roi_i) = sum(ROI_idx==roi_i);
end
xlim([time_vec(1) time_vec(end)]); ylim(ylimit);
xlabel('Time (ms)'); ylabel(settings.ylab); title(settings.title);
set(gca,'FontSize',14,'Box','off')

% add inset with the number of electrodes per region
inset_pos = get(gca,'Position'); 
inset_pos = inset_pos.*[1 1 0.055+0.03*n_reg 0.2] + inset_pos([3 4 3 4]).*[0.92-0.03*n_reg 0.7 0 0];
axes('Units','Normalized','Position',inset_pos);
b = bar(e_nums, 0.7,'EdgeColor','none','FaceAlpha',0.8,'FaceColor','flat'); b.CData = settings.reg_cmap;
ylim_cur = ylim;
text(1:n_reg,repmat(1.15*ylim_cur(2),1,n_reg),string(e_nums),'FontSize',10,'HorizontalAlignment','center');
xlim([0.5 n_reg+0.5]);ylabel('# elec');set(gca,'XTicklabels',settings.ROInames,'FontSize',11);box off;
end

function dval = cohens_d(x1,x2) % for independent samples
n1 = length(x1); n2 = length(x2);
pooledSD =  sqrt(((n1-1)*var(x1) + (n2-1)*var(x2)) / (n1 + n2 - 2)); % pooled Standard Deviation
dval =  (mean(x1) - mean(x2)) / pooledSD; % Cohen's d
end