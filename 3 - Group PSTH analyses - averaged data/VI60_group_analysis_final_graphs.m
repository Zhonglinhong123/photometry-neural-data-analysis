%% Create group-level PSTH for neural activity in VI60 lever press task

% ** What you can expect this script to do:
%  Generate group average graphs of photometry responses (rewarded and unrewarded port entry, lever press) as PSTH

% ** What you need to customize:
%       Import/export data filenames
%       Axis limits for figures (can't be larger than original analysis window)

% ** Assumptions:
%   Raw data have been previously extracted from TDT data files and photometry signal has been fitted and z-score normalized using TDT_photometry_data_extraction_script
%   Data has been previously analyzed to generate subject-level PSTH for photometry with VI60_lever_batch_analysis script

% ** Requires one additional function: 
%       opacity (to generate fills for shaded error bars for graphs)

%% Prepare workspace

clear all
close all

%% Specify input and output filenames and figure names

input_datafile_name = 'VI60_group_analyzed_data';
figname1 = 'VI60_Rew_v_NoRew_PortEntry_PSTH_group_data';
figname2 = 'VI60_active_lever_PSTH_group_data';

%% Specify axis limits for figures 

% Figure 1 - photometry data, PSTH for rewarded/unreward entries 
xmin1 = 10;      % time in seconds - can't be larger than original analysis window
xmax1 = 10;
ymin1 = -1;      % Z-score
ymax1 = 3;

% Figure 2 - photometry data, PSTH for active lever
xmin2 = 1;      % time in seconds - can't be larger than original analysis window
xmax2 = 1;
ymin2 = -1;      % Z-score
ymax2 = 3;

%% Load previously analyzed individual data

load (input_datafile_name)   
                                              
% Key variables for import: 
%       Photometry_PSTH - Structure array with individual-level data in 9 fields:
%       (1) subjectID (2) original datafile (3) sampling rate
%       (4) Rewarded Entry PSTH and (5) SEM
%       (6) Unrewarded Entry PSTH and (7) SEM
%       (8) Active lever PSTH and (9) SEM

%       sec_preEvent and sec_postEvent - time window of original analysis
        originalSecPrev = sec_preEvent;
        originalSecPost = sec_postEvent; 
        
%       Number of subjects used for group analysis
        nSubjects = size(Photometry_PSTH, 2);

% Report an error if x-axis limits exceed time window for original analysis
if xmin1 > originalSecPrev | xmin2> originalSecPrev
    error ('The x-axis minimum limits for group graphs exceed the original analysis window for individual data.');
end

if xmax1 > originalSecPost | xmax2 > originalSecPost
    error ('The x-axis maximum limits for group graphs exceed the original analysis window for individual data.');
end


%% Average individual PSTHs into a group mean vector for plotting - Reward/No reward - Photometry

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.RewardEntry_PSTH}';
group_photometry_RewardEntry_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.UnrewardEntry_PSTH}';
group_photometry_UnrewardEntry_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% make PSTH for rewarded licks
Group_err_Rew = (nanstd(group_photometry_RewardEntry_PSTH))/sqrt(size(group_photometry_RewardEntry_PSTH, 1));
GroupPSTH_Rew = nanmean(group_photometry_RewardEntry_PSTH);

% make PSTH for unrewarded licks
Group_err_Unrew = (nanstd(group_photometry_UnrewardEntry_PSTH))/sqrt(size(group_photometry_UnrewardEntry_PSTH, 1));
GroupPSTH_Unrew = nanmean(group_photometry_UnrewardEntry_PSTH);


%% Make final graph - Reward/No reward - Photometry

% Create graph
Figure1 = figure;
hold on

nTsPrev = round (xmin1 * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (xmax1 * samplingRate);
totalTs = nTsPrev + nTsPost;                     % sets time axis for upcoming graph based on these values
increment = (xmin1 + xmax1) / totalTs;
timeAxis = (-1 * xmin1) : increment : xmax1;
ntimeAxis = size (timeAxis, 2);
axis_start = round ((originalSecPrev-xmin1)* samplingRate);
if axis_start == 0
    axis_start = 1;
end
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
ylabel('Z-score', 'FontSize', 20);
legend ([h1 h2], 'Reward', 'No reward', 'orientation', 'vertical', 'Location', 'NorthWest');
legend BOXOFF;
title('VI60 - Rewarded and Unrewarded Port Entries');


%% Average individual PSTHs into a group mean vector for plotting - Active lever - Photometry

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.Active_Lever_PSTH}';
group_photometry_Active_Lever_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% make PSTH for early rewarded licks
Group_err_Lever = (nanstd(group_photometry_Active_Lever_PSTH))/sqrt(size(group_photometry_Active_Lever_PSTH, 1));
GroupPSTH_Lever = nanmean(group_photometry_Active_Lever_PSTH);


%% Make final graph - Active lever press

% Create graph
Figure2 = figure;
hold on
   
nTsPrev = round (xmin2 * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (xmax2 * samplingRate);
totalTs = nTsPrev + nTsPost;                          % sets time axis for upcoming graph based on these values
increment = (xmin2 + xmax2) / totalTs;
timeAxis = (-1 * xmin2) : increment : xmax2;
ntimeAxis = size (timeAxis, 2);
axis_start = round ((originalSecPrev-xmin2)* samplingRate);
if axis_start == 0
    axis_start = 1;
end
axis_end = axis_start + ntimeAxis - 1;

Err1Pos = GroupPSTH_Lever + Group_err_Lever;
Err1Neg = GroupPSTH_Lever - Group_err_Lever;

o = 0.5;

r1 = 24;
g1 = 86;
b1 = 9;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos(1, axis_start:axis_end), fliplr(Err1Neg(1, axis_start:axis_end))], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,GroupPSTH_Lever(1, axis_start:axis_end),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

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
legend ([h1], 'Active Lever', 'orientation', 'vertical', 'Location', 'NorthWest');
legend BOXOFF;
title('VI60 - Active Lever');

%% Save graphs

saveas (Figure1, figname1);
saveas (Figure2, figname2);

clearvars
    