# PersistentViewing_paper
Matlab code for analysis and visualization for: Vishne, Gerber, Knight and Deouell, Cell Reports 2023 “Distinct Ventral Stream and Prefrontal Cortex Representational Dynamics during Sustained Conscious Visual Perception” (biorxiv DOI, to be updated when formally published): https://doi.org/10.1101/2022.08.02.502469

Please cite if you use this! 🙏

## Code structure:

The paper has five sections, using four distinct analysis approaches, accordingly, the code is divided into four folders:
1. Univariate analyses - results section 1: Figure 1 and Supplementary Figure 1
2. Multivariate state-space analyses - results section 2: Figure 2 and Supplementary Figures 2-3
3. Decoding category information - containing two parts:
 1. Results section 3, focusing on long stimuli durations (more than 900ms): Figure 3 and Supplementary Figures 4-6
 2. Results section 4, comparing decoding in different durations: Figure 4 and Supplementary Figure 7
4. Single-exemplar analyses - results section 5: Figure 5 and Supplementary Figures 8-11.

For each section the code is divided into a X_calc.m script and X_fig.m script. The first is used to run all analyses and save the results, and the second loads this information and plots the Figures in the manuscript. In some folders there are additional functions used for visualization or analysis specific to those sections.

The analysis scripts (X_calc.m) call data produced with the script specific_segs.m. This script uses basic_segs_1000Hz.mat or basic_segs_200Hz.mat (see details below), to create the specific segment structure required for each analysis (merging the segments to one ‘mega-patient’, explanation of the relevant settings is inside the function). Basic_segs.mat and other data can be found on: OSF (repository name: PersistentViewingECoG, DOI 10.17605/OSF.IO/4HXPW).

The figures can all be reproduced by using this order: specific_segs.m -> X_calc.m -> X_fig.m, but I also uploaded to OSF the intermediate output from X_calc.m since permutation results may be slightly different on your computer.


## Additional code:

The “heavy lifting” for the analyses in results sections 3-5 is done using toolboxes I built for time-resolved analyses (links below). All three repositories include a lot of functionality that wasn’t used in the paper, so you are welcome to explore:
- Decoding: **iEEG_decoding_minitoolbox** (https://github.com/NeuroGal/iEEG_decoding_minitoolbox) - this toolbox provides a wrapper function which calls MVPA_Light (Treder, MS Frontiers in Neuroscience (2020) https://doi.org/10.3389/fnins.2020.00289) to do the decoding, but the statistics are all implemented separately in a way specific to intracranial EEG (iEEG\ECoG) or spike data when you have one 'super subject'. The decoding_stats.m function calls different functions implemented in the time_resolved_stats below.
- Representational similarity analysis: **Gals_RSA_toolbox** (https://github.com/NeuroGal/Gals_RSA_toolbox) - this toolbox includes three parts: 1. Code to generate the Representational Dissimilarity Matrix (RDM), including code for noise normalization prior to the calculation (also includes plotting). 2. Code for exemplar discriminability index (also includes plotting). 3. Code for calculating single-exemplar information reliability across stimulus repetitions \ time. For this manuscript I use only parts 1 & 3. Parts 2 & 3 include statistics which use the functions implemented in time_resolved_stats below.
- Statistical functions (used by both): **time_resolved_stats** (https://github.com/NeuroGal/time_resolved_stats) - implements cluster-based permutations, max-statistic control for multiple comparisons and more.

Visualization is done using custom code in each folder and general functions which are used in multiple sections. These are included in a folder marked as “Helper functions” (which also includes a bit of non-plotting related functions). The folder includes also:
- Tight_subplot.m written by Pekka Kumpulainen (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
- Varplot.m written by Edden Gerber (https://www.mathworks.com/matlabcentral/fileexchange/71632-varplot)
In both cases I put the version I worked with in the folder.

Dimensionality reduction and visualization of the multivariate response (Figures 2A, Supplementary Figures 2A and 11A) is done using my repository **state_space_plot** (https://github.com/NeuroGal/state_space_plot). This implements both plotting of each time-point separately (as in Supplementary Figure 11A), and plotting all together (state-space trajectories, as in Figures 2A and Supplementary Figures 2A). There are also a ton of additional options I didn’t use in any Figure, so I encourage you to check it out yourself.


## Data:
Deposited on OSF: repository name - PersistentViewingECoG, DOI - 10.17605/OSF.IO/4HXPW
Basic_segs_200Hz.mat and Basic_segs_1000Hz.mat each contain the following variables:

- Segs_basic: cell(10,1) with each cell containing an array of (stimulus x electrode x time) of HFA activity of each patient, after downsampling, followed by moving average smoothing with 50ms and baseline correction by subtracting 300ms prior to stimulus onset. The results are very similar if you do the downsampling at a later stage, but since this was the order it was performed in the manuscript I enclosed both a 200Hz and a 1000Hz version.
- Trial_info_basic: cell(10,1) with the information corresponding to the stimuli dimension in segs - stimulus category, id, duration and ISI.
- Time_vec_basic: time (around stimulus onset) corresponding to the time dimension of segs (changes when you change the sampling rate).
- Cat_names: names of the categories (trial_info_basic has numbers)
- Elec_info: table with electrode properties (including all patients) - fields:
  - Patient (numbers 1-10)
  - IsResponsive (only responsive electrodes were used for more analyses)
  - RespSign (+\-, increase or decrease from baseline)
  - RespPerCat_FWOA (responsive per category, according to the order faces,watches,objects,animals)
  - IsSelective (category selectivity)
  - ROI (Occ,VT,Par,PFC,SM,LT, see manuscript for the acronyms)
  - ROISplit (same as ROI but separating Occ to retinotopic\non-retinotopic and PFC to OFC and LPFC)
  - ROIDetail (more detailed anatomical location)
  - Hemisphere
  - MNIcoord (x,y,z coordinates on MNI152 template from freesurfer, this was used to image the brains in the manuscript figures).

Settings needed for the plotting are in: figure_settings.mat (mostly includes color maps).

As written above - this is enough to get you to all figures (with one exception noted below) using the following order: specific_segs.m -> X_calc.m -> X_fig.m, but the intermediate output from the different analysis scripts (X_calc.m of each section) is also found on OSF. 

For 2-repetition data used in the single-exemplar analyses specific_segs.m contains a randomization step. If you want to use the exact segments I used in the manuscript that can be found on OSF with the name paper50ms_XXX.mat or paper100ms_XXX.mat.

For Supplementary Figure 6e, comparing decoding with\without saccades in the relevant window use segs_ofc_saccade_control.mat, which contains saccade information not found in basic_segs.mat. This is used in the script ofc_saccade_control.m (in the decoding analyses folder) which has more explanation about this analysis.


## A few final details:
- For imaging brains (Figures 1C, Supplementary Figures 1D,F, 5A-C, 6A-C, 7A and 8A) I use the coordinates in ‘elec_info’ and custom scripts - I didn’t have time to sort them out, but you are welcome to contact me for this. The scripts are based on Edden Gerber’s wonderful repository: https://edden-gerber.github.io/vis_toolbox/ 
- HFA data was computed by filtering eight 10 Hz bands from 70-150 Hz, using pop_eegfiltnew.m from EEGLab then taking the abs(hilbert()) of each output. Each band was normalized by dividing all time points by the mean amplitude throughout the recording and then the signals were averaged (see the STAR methods for more detail).
- Electrode responsiveness & category selectivity was computed using 200 ms windows with unsmoothed data, and are therefore included as part of ‘elec_info’ (in basic_segs.mat). The details are explained in STAR Methods of the manuscript (conclusions hold with other definitions too).
- Basic_segs_1000Hz is in single format due to size limitations on OSF. The calculations were done on double data. I assume it wouldn’t matter much. If needed the interim calculations are also on OSF.


Gal Vishne, June 2023

Twitter: @neuro_gal
