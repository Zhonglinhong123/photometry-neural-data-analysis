%% Create group-level PSTH for behavior and neural activity in a lickometer task

% ** What you can expect this script to do:
%  Generate group average graphs of behavioral data and photometry responses as PSTH

% ** What you need to customize:
%       Import/export data filenames

% ** Assumptions:
%   Raw data have been previously extracted from TDT data files and photometry signal has been fitted and z-score normalized using TDT_photometry_data_extraction_script
%   Data has been previously analyzed to generate subject-level PSTH for behavior and photometry with Lickometer_batch_analysis_github script

% ** Requires one additional function: 
%       opacity (to generate fills for shaded error bars for graphs)

%% Prepare workspace

clear all
close all

%% Specify input and output filenames and figure names

input_datafile_name = 'Lickometer_group_analyzed_data';
figname1 = 'Lickometer_behavior_lick_rate_PSTH_group_data';
figname2 = 'Lickometer_photometry_rewarded-unrewarded_lick_PSTH_group_data';
figname3 = 'Lickometer_photometry_rewarded_licks_early-late_group_data';
figname4 = 'Lickometer_photometry_unrewarded_licks_early-late_group_data';

%% Specify axis limits for figures 

% Figure 1 - behavioral data, lick rate for rewarded and unrewarded licks
xmin1 = 2;      % time in seconds - can't be larger than original analysis window
xmax1 = 4;
ymin1 = 0;      % Licks/sec
ymax1 = 15;


% Figure 2 - photometry data, PSTH for rewarded/unreward licks (all trials)
xmin2 = 2;      % time in seconds - can't be larger than original analysis window
xmax2 = 4;
ymin2 = -1;      % Z-score
ymax2 = 3;

% Figure 3 - photometry data, PSTH for early/late rewarded trials
xmin3 = 2;      % time in seconds - can't be larger than original analysis window
xmax3 = 4;
ymin3 = -1;      % Z-score
ymax3 = 6;

% Figure 4 - photometry data, PSTH for early/late unrewarded trials
xmin4 = 2;      % time in seconds - can't be larger than original analysis window
xmax4 = 4;
ymin4 = -1;      % Z-score
ymax4 = 6;

%% Load previously analyzed individual data

load (input_datafile_name)   
                                              
% Key variables for import: 
%       Photometry_PSTH - Structure array with 15 fields:
%       (1) subjectID (2) original datafile (3) sampling rate
%       (4) Rewarded Lick PSTH and (5) SEM
%       (6) Unrewarded Lick PSTH and (7) SEM
%       (8) Early rewarded lick PSTH and (9) SEM
%       (10) Late rewarded lick PSTH and (11) SEM
%       (12) Early unrewarded lick PSTH and (13) SEM
%       (14) Late unrewarded lick PSTH and (15) SEM

%       Behavior PSTH - Structure array with 6 fields:
%       (1) subjectID (2) original datafile
%       (3) Rewarded lick PSTH and (4) SEM
%       (5) Unrewarded lick PSTH and (6) SEM

%       sec_preEvent and sec_postEvent - time window of original analysis
        originalSecPrev = sec_preEvent;
        originalSecPost = sec_postEvent; 
        
%       Number of subjects used for group analysis
nSubjects = size(Behavior_PSTH, 2);

% Report an error if x-axis limits exceed time window for original analysis
if xmin1 > originalSecPrev | xmin2> originalSecPrev | xmin3> originalSecPrev | xmin4 > originalSecPrev
    error ('The x-axis minimum limits for group graphs exceed the original analysis window for individual data.');
end

if xmax1 > originalSecPost | xmax2 > originalSecPost | xmax3 > originalSecPost | xmax4 > originalSecPost
    error ('The x-axis maximum limits for group graphs exceed the original analysis window for individual data.');
end

%% Average individual PSTHs into a group mean vector for plotting - Reward/No reward - Behavior

% Extract individual PSTH vectors from structure array
tempCell = {Behavior_PSTH.rewLick_PSTH}';
group_behavior_rewLick_PSTH = cell2mat(tempCell);

tempCell = {Behavior_PSTH.UnrewLick_PSTH}';
group_behavior_UnrewLick_PSTH = cell2mat(tempCell);

% make group PSTH for rewarded licks
Group_err_Rew_beh = (nanstd(group_behavior_rewLick_PSTH))/sqrt(size(group_behavior_rewLick_PSTH, 1));
GroupPSTH_Rew_beh = nanmean(group_behavior_rewLick_PSTH);

% make group PSTH for unrewarded licks
Group_err_Unrew_beh = (nanstd(group_behavior_UnrewLick_PSTH))/sqrt(size(group_behavior_UnrewLick_PSTH, 1));
GroupPSTH_Unrew_beh = nanmean(group_behavior_UnrewLick_PSTH);

%% Make final graph - Reward/No reward - Behavior 

% Create graph
Figure1 = figure;
hold on

lickTotalSec = xmin1+xmax1;
beh_timeAxis = (-1 * xmin1) : beh_time_resolution : xmax1;
nBehtimeAxis = size (beh_timeAxis, 2);
axis_start = ((originalSecPrev - xmin1)/beh_time_resolution)+1;
axis_end = axis_start + nBehtimeAxis - 1;

Err1Pos = GroupPSTH_Rew_beh + Group_err_Rew_beh;
Err1Neg = GroupPSTH_Rew_beh - Group_err_Rew_beh;
Err2Pos = GroupPSTH_Unrew_beh + Group_err_Unrew_beh;
Err2Neg = GroupPSTH_Unrew_beh - Group_err_Unrew_beh;

o = 0.5;

r2 = 175;
g2 = 175;
b2 = 175;
rgb2_o = opacity (o, r2, g2, b2);

fill([beh_timeAxis, fliplr(beh_timeAxis)],[Err2Pos(1, axis_start:axis_end), fliplr(Err2Neg(1, axis_start:axis_end))], rgb2_o, 'EdgeColor', 'none');
h2= plot (beh_timeAxis,GroupPSTH_Unrew_beh(1, axis_start:axis_end),'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

r1 = 0;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);

fill([beh_timeAxis, fliplr(beh_timeAxis)],[Err1Pos(1, axis_start:axis_end), fliplr(Err1Neg(1, axis_start:axis_end))], rgb1_o, 'EdgeColor', 'none');
h1= plot (beh_timeAxis,GroupPSTH_Rew_beh(1, axis_start:axis_end),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);


%labels, legend, make pretty, size
line('XData', [(-1 * xmin1),xmax1], 'YData', [0 0], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin1, ymax1], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
xlim ([(-1 * xmin1),xmax1]);
ylim ([ymin1, ymax1]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 18, ...
  'LineWidth'   , 2         );
txt = ['n = ' num2str(nSubjects)];
text_xpos = xlim*0.7;   % adjust scaling factor to adjust text position
text_ypos = ylim*0.9;   
text(text_xpos(2), text_ypos(2), txt, 'FontSize',18);

xlabel('Time (s)', 'FontSize', 20);
ylabel('Licks/s', 'FontSize', 20);
legend ([h1 h2], 'Reward', 'No reward', 'orientation', 'vertical', 'Location', 'NorthWest');
legend BOXOFF;
title('Behavior - Lick rate');


%% Average individual PSTHs into a group mean vector for plotting - Reward/No reward - Photometry

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.RewardLick_PSTH}';
group_photometry_RewardLick_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.UnrewardLick_PSTH}';
group_photometry_UnrewardLick_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% make PSTH for rewarded licks
Group_err_Rew = (nanstd(group_photometry_RewardLick_PSTH))/sqrt(size(group_photometry_RewardLick_PSTH, 1));
GroupPSTH_Rew = nanmean(group_photometry_RewardLick_PSTH);

% make PSTH for unrewarded licks
Group_err_Unrew = (nanstd(group_photometry_UnrewardLick_PSTH))/sqrt(size(group_photometry_UnrewardLick_PSTH, 1));
GroupPSTH_Unrew = nanmean(group_photometry_UnrewardLick_PSTH);


%% Make final graph - Reward/No reward - Photometry

% Create graph
Figure2 = figure;
hold on

nTsPrev = round (xmin2 * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (xmax2 * samplingRate);
totalTs = nTsPrev + nTsPost;                     % sets time axis for upcoming graph based on these values
increment = (xmin2 + xmax2) / totalTs;
timeAxis = (-1 * xmin2) : increment : xmax2;
ntimeAxis = size (timeAxis, 2);
axis_start = round ((originalSecPrev-xmin2)* samplingRate);
axis_end = axis_start + ntimeAxis - 1;

Err1Pos = GroupPSTH_Rew + Group_err_Rew;
Err1Neg = GroupPSTH_Rew - Group_err_Rew;
Err2Pos = GroupPSTH_Unrew + Group_err_Unrew;
Err2Neg = GroupPSTH_Unrew - Group_err_Unrew;

o = 0.5;

r2 = 75;
g2 = 75;
b2 = 75;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[Err2Pos(1, axis_start:axis_end), fliplr(Err2Neg(1, axis_start:axis_end))], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,GroupPSTH_Unrew(1, axis_start:axis_end),'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

r1 = 0;
g1 = 0;
b1 = 255;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos(1, axis_start:axis_end), fliplr(Err1Neg(1, axis_start:axis_end))], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,GroupPSTH_Rew(1, axis_start:axis_end),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
line('XData', [(-1 * xmin2),xmax2], 'YData', [0 0], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin2, ymax2], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
xlim ([(-1 * xmin2),xmax2]);
ylim ([ymin2, ymax2]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 18, ...
  'LineWidth'   , 2         );
txt = ['n = ' num2str(nSubjects)];
text_xpos = xlim*0.7;   % adjust scaling factor to adjust text position
text_ypos = ylim*0.9;   
text(text_xpos(2), text_ypos(2), txt, 'FontSize',18);

xlabel('Time (s)', 'FontSize', 20);
ylabel('Z-score', 'FontSize', 20);
legend ([h1 h2], 'Reward', 'No reward', 'orientation', 'vertical', 'Location', 'NorthWest');
legend BOXOFF;
title('Photometry data, all trials');


%% Average individual PSTHs into a group mean vector for plotting - Reward early/late - Photometry

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.early_RewLick_PSTH}';
group_photometry_early_RewLick_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.late_RewLick_PSTH}';
group_photometry_late_RewLick_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% make PSTH for early rewarded licks
Group_err_RE = (nanstd(group_photometry_early_RewLick_PSTH))/sqrt(size(group_photometry_early_RewLick_PSTH, 1));
GroupPSTH_RE = nanmean(group_photometry_early_RewLick_PSTH);

% make PSTH for late rewarded licks
Group_err_RL = (nanstd(group_photometry_late_RewLick_PSTH))/sqrt(size(group_photometry_late_RewLick_PSTH, 1));
GroupPSTH_RL = nanmean(group_photometry_late_RewLick_PSTH);

%% Make final graph - Reward early/late - Photometry

% Create graph
Figure3 = figure;
hold on
   
nTsPrev = round (xmin3 * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (xmax3 * samplingRate);
totalTs = nTsPrev + nTsPost;                          % sets time axis for upcoming graph based on these values
increment = (xmin3 + xmax3) / totalTs;
timeAxis = (-1 * xmin3) : increment : xmax3;
ntimeAxis = size (timeAxis, 2);
axis_start = round ((originalSecPrev-xmin3)* samplingRate);
axis_end = axis_start + ntimeAxis - 1;

Err1Pos = GroupPSTH_RE + Group_err_RE;
Err1Neg = GroupPSTH_RE - Group_err_RE;
Err2Pos = GroupPSTH_RL + Group_err_RL;
Err2Neg = GroupPSTH_RL - Group_err_RL;

o = 0.75;

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[Err2Pos(1, axis_start:axis_end), fliplr(Err2Neg(1, axis_start:axis_end))], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,GroupPSTH_RL(1, axis_start:axis_end),'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3, 'LineStyle', ':');

o = 0.5;

r1 = 0;
g1 = 0;
b1 = 255;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos(1, axis_start:axis_end), fliplr(Err1Neg(1, axis_start:axis_end))], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,GroupPSTH_RE(1, axis_start:axis_end),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
line('XData', [(-1 * xmin3),xmax3], 'YData', [0 0], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin3, ymax3], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
xlim ([(-1 * xmin3),xmax3]);
ylim ([ymin3, ymax3]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 18, ...
  'LineWidth'   , 2         );
txt = ['n = ' num2str(nSubjects)];
text_xpos = xlim*0.7;   % adjust scaling factor to adjust text position
text_ypos = ylim*0.9;   
text(text_xpos(2), text_ypos(2), txt, 'FontSize',18);

xlabel('Time (s)', 'FontSize', 20);
ylabel('Z-score', 'FontSize', 20);
legend ([h1 h2], 'Reward (early)', 'Reward (late)', 'orientation', 'vertical', 'Location', 'NorthWest');
legend BOXOFF;
title('Photometry data, early vs late rewarded trials');


%% Average individual PSTHs into a group mean vector for plotting - No reward early/late - Photometry

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.early_UnrewLick_PSTH}';
group_photometry_early_UnrewLick_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.late_UnrewLick_PSTH}';
group_photometry_late_UnrewLick_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;


% make PSTH for early unrewarded licks
Group_err_UE = (nanstd(group_photometry_early_UnrewLick_PSTH))/sqrt(size(group_photometry_early_UnrewLick_PSTH, 1));
GroupPSTH_UE = nanmean(group_photometry_early_UnrewLick_PSTH);

% make PSTH for late unrewarded licks
Group_err_UL = (nanstd(group_photometry_late_UnrewLick_PSTH))/sqrt(size(group_photometry_late_UnrewLick_PSTH, 1));
GroupPSTH_UL = nanmean(group_photometry_late_UnrewLick_PSTH);

%% Make final graph - No reward early/late - Photometry 

% Create graph
Figure4 = figure;
hold on

nTsPrev = round (xmin4 * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (xmax4 * samplingRate);
totalTs = nTsPrev + nTsPost;                     % sets time axis for upcoming graph based on these values
increment = (xmin4 + xmax4) / totalTs;
timeAxis = (-1 * xmin4) : increment : xmax4;
ntimeAxis = size (timeAxis, 2);
axis_start = round ((originalSecPrev-xmin4)* samplingRate);
axis_end = axis_start + ntimeAxis - 1;

Err1Pos = GroupPSTH_UE + Group_err_UE;
Err1Neg = GroupPSTH_UE - Group_err_UE;
Err2Pos = GroupPSTH_UL + Group_err_UL;
Err2Neg = GroupPSTH_UL - Group_err_UL;

o = 0.5;

r2 = 150;
g2 = 150;
b2 = 150;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[Err2Pos(1, axis_start:axis_end), fliplr(Err2Neg(1, axis_start:axis_end))], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,GroupPSTH_UL(1, axis_start:axis_end),'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3, 'LineStyle', ':');

r1 = 75;
g1 = 75;
b1 = 75;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos(1, axis_start:axis_end), fliplr(Err1Neg(1, axis_start:axis_end))], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,GroupPSTH_UE(1, axis_start:axis_end),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
line('XData', [(-1 * xmin4),xmax4], 'YData', [0 0], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin4, ymax4], 'LineStyle', '--', 'LineWidth', 1, 'Color','k')
xlim ([(-1 * xmin4),xmax4]);
ylim ([ymin4, ymax4]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 18, ...
  'LineWidth'   , 2         );
txt = ['n = ' num2str(nSubjects)];
text_xpos = xlim*0.7;   % adjust scaling factor to adjust text position
text_ypos = ylim*0.9;   
text(text_xpos(2), text_ypos(2), txt, 'FontSize',18);

xlabel('Time (s)', 'FontSize', 20);
ylabel('Z-score', 'FontSize', 20);
legend ([h1 h2], 'No reward (early)', 'No reward (late)', 'orientation', 'vertical', 'Location', 'NorthWest');
legend BOXOFF;
title('Photometry data, early vs late unrewarded trials');


%% Save graphs

saveas (Figure1, figname1);
saveas (Figure2, figname2);
saveas (Figure3, figname3);
saveas (Figure4, figname4);

clearvars
    