%% Create subject-level PSTH for behavior and neural activity in a lickometer task where sucrose is intermittently available

% ** What you can expect this script to do:
%  Select timestamps for rewarded and unrewarded licks that are separated in time to avoid re-sampling the same data twice
%  Generate batch graphs of behavioral data and photometry responses as PSTH for QC and exploratory visualization of individual variability
%  Key output data stored in three structure arrays: 
 
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

%       Behavior_ts - Structure array with 6 fields:
%       (1) subjectID (2) original datafile
%       (3) Lick timestamps (4) Reward timestamps
%       (5) Rewarded lick timestamps (6) Unrewarded lick timestamps (processed subset)


% ** What you need to customize:
%       Import/export data filenames
%       Source of behavioral data - separate matlab file or from TTL inputs through TDT system (will default to TTL inputs through TDT system if behavior filenames are not supplied)
%       Identify which TTL inputs correspond to licks and rewards

% ** Assumptions:
%   Raw data have been previously extracted from TDT data files and photometry signal has been fitted and z-score normalized using TDT_photometry_data_extraction_script
%   See below for expected input data format

% ** Requires two additional functions: 
%       processPhotDataRow_normDat (to align photometry data across trials to generate PSTH, credit: Tom Davidson)
%       opacity (to generate fills for shaded error bars for graphs)

%% Prepare workspace

clear all
close all

%% Specify input and output filenames and figure names and number of rows/columns for graphs with individual subplots for each subject, 

input_datafile_name = 'Lickometer_group_raw_data';
output_datafile_name = 'Lickometer_group_analyzed_data';
figname1 = 'Lickometer_behavior_lick_rate_PSTH_individual_data';
figname2 = 'Lickometer_photometry_rewarded-unrewarded_lick_PSTH_individual_data';
figname3 = 'Lickometer_photometry_rewarded_licks_early-late_individual_data';
figname4 = 'Lickometer_photometry_unrewarded_licks_early-late_individual_data';
subject_plot_rows = 1;                                                      
subject_plot_cols = 3; 

% Optional - specify separate files with behavioral data
behavior_files = {'FP25_Lickometer_030217.mat'; 'FP41_Lickometer_030217.mat'; 'FP42_Lickometer_030217.mat'};

% Specify PSTH time window
sec_preEvent = 5;
sec_postEvent = 5;

%% Load previously extracted and normalized photometry data

load (input_datafile_name)                                         
% Expected variables for import: Photometry_zScore structure array with >= 6 fields
% Field 1: subject IDs
% Field 2: original datafile names
% Field 3: sampling rate
% Field 4: fitted, z-score normalized calcium-dependent photometry signal
% Field 5: timestamps for photometry signal
% Field 6 and higher: TTL timestamps corresponding to external events

%% Run analysis loop for each subject

for c = 1:length(subject_IDs)
    
%% Load behavioral data - code for option 1 executes if behavioral data file names are provided above; otherwise defaults to option 2

% Option 1 - use this code if there are separate matlab file(s) with behavioral data
if numel(behavior_files) >=1 
    behfile = char(behavior_files(c));
    load (behfile);
    timezero_ts = Photometry_zScore(c).TTL_event1;
    Rewards_ts = reward + timezero_ts;
    Licks_ts = licks + timezero_ts;
else
% Option 2 - use this code if behavioral events are captured in the same TDT datafile as photometry signal
    Rewards_ts = Photometry_zScore(c).TTL_event1;
    Licks_ts = Photometry_zScore(c).TTL_event2;
end

%% Load photometry data

zDat = Photometry_zScore(c).photometry_zScore_data;
samplingRate = Photometry_zScore(c).photometry_samplingrate;

%% Find first lick following reward

nReward_ts = size(Rewards_ts, 1);
nLicks = size(Licks_ts, 1);
RewLicks_ts = NaN(nReward_ts, 1);
for iRew = 1:nReward_ts;
    thisRew = Rewards_ts(iRew);
    for iLicks = 1:nLicks;
        thisLick = find(Licks_ts>=thisRew, 1, 'first');
        RewLicks_ts(iRew,1) = Licks_ts(thisLick,1);      % rewarded lick = 1st lick after reward
    end
end


%% Find unrewarded licks that don't overlap with rewarded lick PSTH or other unrewarded lick PSTH to avoid resampling data

PreRewCriterion = 5;        % Sec before reward to exclude                                     
PostRewCriterion = 5;       % Sec after reward to exclude

UnrewLicks_ts = Licks_ts;
nUnrewLicks = size(UnrewLicks_ts, 1);
nRewLicks = size(RewLicks_ts, 1);

for iRew = 1:nRewLicks;                                                     % Must fall outside of criterion times before and after rewarded licks
    thisRewLick = RewLicks_ts(iRew, 1);
    for iUnRew = 1:nUnrewLicks;
        if UnrewLicks_ts(iUnRew, 1) > (thisRewLick - PreRewCriterion) & UnrewLicks_ts(iUnRew, 1) < (thisRewLick + PostRewCriterion);
           UnrewLicks_ts(iUnRew, 1) = NaN;
        end
    end
end
UnrewLicks_ts = UnrewLicks_ts(all(~isnan(UnrewLicks_ts), 2),:);             % Unrewarded PE times - BUT these may occur in bursts which are problematic for PSTH averaging


PostUnrewCriterion = 5;                                                     % Remove unrewarded licks that occur within criterion sec of prev unrewarded PE (to avoid duplicating same neural data with slight time shift in PSTH)
nUnrewLicks = size(UnrewLicks_ts, 1);
for iUnRew = 1:nUnrewLicks;
    remainder = nUnrewLicks - iUnRew;
    for i = 1:remainder;
        if UnrewLicks_ts(iUnRew+i, 1) <= (UnrewLicks_ts(iUnRew, 1) + PostUnrewCriterion);
           UnrewLicks_ts(iUnRew+i, 1) = NaN;
        end
    end
end

UnrewLicks_ts = UnrewLicks_ts(all(~isnan(UnrewLicks_ts), 2),:);             % This is processed set of lick times that meet all criteria for being unrewarded and are separated in time


%% Calculate lick rate histograms (behavior readout)

lickTotalSec = sec_preEvent+sec_postEvent;
beh_time_resolution = 0.2;
nBins = lickTotalSec/beh_time_resolution+1;
nLicks = size(Licks_ts, 1);

% Rewarded lick histogram
nRewLicks = size(RewLicks_ts, 1);
Rew_lick_array = zeros(nRewLicks, nBins);
for iRew = 1:nRewLicks;
    trialStart = RewLicks_ts(iRew, 1)-sec_preEvent;
    for iBin = 1:nBins;
        thisBin = Rew_lick_array(1, iBin);
        binStart = trialStart+((iBin-1)*beh_time_resolution);
        binEnd = trialStart+(iBin*beh_time_resolution);      
        for iLick = 1:nLicks;
            thisLick = Licks_ts(iLick,1);
            if thisLick >= binStart & thisLick <binEnd;
                Rew_lick_array (iRew, iBin) = Rew_lick_array (iRew, iBin)+1;
            end
        end
    end
end

beh_rewLick_PSTH = (mean(Rew_lick_array))/beh_time_resolution;
beh_rewLick_err = (std(Rew_lick_array)/sqrt(size(Rew_lick_array,1)))/beh_time_resolution;


% Unrewarded lick histogram
nUnrewLicks = size(UnrewLicks_ts, 1);
Unrew_lick_array = zeros(nUnrewLicks, nBins);
for iUnrew = 1:nUnrewLicks;
    trialStart = UnrewLicks_ts(iUnrew, 1)-sec_preEvent;
    for iBin = 1:nBins;
        thisBin = Unrew_lick_array(1, iBin);
        binStart = trialStart+((iBin-1)*beh_time_resolution);
        binEnd = trialStart+(iBin*beh_time_resolution);      
        for iLick = 1:nLicks;
            thisLick = Licks_ts(iLick,1);
            if thisLick >= binStart & thisLick <binEnd;
                Unrew_lick_array (iUnrew, iBin) = Unrew_lick_array (iUnrew, iBin)+1;
            end
        end
    end
end

beh_UnrewLick_PSTH = (mean(Unrew_lick_array))/beh_time_resolution;
beh_UnrewLick_err = (std(Unrew_lick_array)/sqrt(size(Unrew_lick_array,1)))/beh_time_resolution;

%% Plot lick rate histograms for each subject (behavior readout)

% Set time window for lick rategraphs
lickTotalSec = sec_preEvent+sec_postEvent;
beh_time_resolution = 0.2;
nBins = lickTotalSec/beh_time_resolution;

beh_timeAxis = (-1 * sec_preEvent) : beh_time_resolution : sec_postEvent;
Err1Pos = beh_rewLick_PSTH + beh_rewLick_err;
Err1Neg = beh_rewLick_PSTH - beh_rewLick_err;
Err2Pos = beh_UnrewLick_PSTH + beh_UnrewLick_err;
Err2Neg = beh_UnrewLick_PSTH - beh_UnrewLick_err;
% Optional - set y-axis limits
% ymin = -1;
% ymax = 1;

if c==1
   batch_licks_beh = figure;
end
figure (batch_licks_beh);
subplot (subject_plot_rows, subject_plot_cols, c);
hold on

% Change RGB values to adjust color, change 'o' value to adjust opacity for shaded error fills
o = 0.5;

r2 = 175;
g2 = 175;
b2 = 175;
rgb2_o = opacity (o, r2, g2, b2);

fill([beh_timeAxis, fliplr(beh_timeAxis)],[Err2Pos, fliplr(Err2Neg)], rgb2_o, 'EdgeColor', 'none');
h2= plot (beh_timeAxis,beh_UnrewLick_PSTH,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

r1 = 0;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);

fill([beh_timeAxis, fliplr(beh_timeAxis)],[Err1Pos, fliplr(Err1Neg)], rgb1_o, 'EdgeColor', 'none');
h1= plot (beh_timeAxis,beh_rewLick_PSTH,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
xlabel('Time (s)', 'FontSize', 10);
ylabel('Licks/s', 'FontSize', 10);
line('XData', [(-1 * sec_preEvent),sec_postEvent], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [0, 10], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
if c == 1
    legend ([h1 h2], 'Reward', 'No reward','orientation', 'vertical', 'Location', 'NorthWest');
    legend BOXOFF;
end
xlim ([(-1 * sec_preEvent),sec_postEvent]);
% ylim ([ymin, ymax]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 10, ...
  'LineWidth'   , 2         );
title(Photometry_zScore(c).subjectID);


%% Make PSTH arrays of photometry data for rewarded and unrewarded licks, then average rows into a mean vector for plotting

% Convert seconds to TDT timestamps
nTsPrev = round (sec_preEvent * samplingRate);
nTsPost = round (sec_postEvent * samplingRate);


% Make PSTH for rewarded licks
nRewLicks = size(RewLicks_ts,1);
PsthArray_Rew = NaN(nRewLicks,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nRewLicks
    thisTime = RewLicks_ts(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_Rew(i,:) = processPhotDataRow_normDat(zDat, thisIndex, nTsPrev, nTsPost);
end
err_Rew = (nanstd(PsthArray_Rew))/sqrt(size(PsthArray_Rew,1));
Psth_Rew = nanmean(PsthArray_Rew);


% Make PSTH for unrewarded licks
nUnrewLicks = size(UnrewLicks_ts,1);
PsthArray_Unrew = NaN(nUnrewLicks,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nUnrewLicks
    thisTime = UnrewLicks_ts(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_Unrew(i,:) = processPhotDataRow_normDat(zDat, thisIndex, nTsPrev, nTsPost);
end
err_Unrew = (nanstd(PsthArray_Unrew))/sqrt(size(PsthArray_Unrew,1));
Psth_Unrew = nanmean(PsthArray_Unrew);


%% Plot rewarded and unrewarded licks - photometry data

totalTs = nTsPrev + nTsPost;
increment = (sec_preEvent + sec_postEvent) / totalTs;
timeAxis = (-1 * sec_preEvent) : increment : sec_postEvent;
Err1Pos = Psth_Rew + err_Rew;
Err1Neg = Psth_Rew - err_Rew;
Err2Pos = Psth_Unrew + err_Unrew;
Err2Neg = Psth_Unrew - err_Unrew;
% Optional - set y-axis limits
% ymin = -1;
% ymax = 1;

if c==1
   batch_licks = figure;
end
figure (batch_licks);
subplot (subject_plot_rows, subject_plot_cols, c);
hold on

% Change RGB values to adjust color, change 'o' value to adjust opacity for shaded error fills
o = 0.5;

r2 = 75;
g2 = 75;
b2 = 75;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[Err2Pos, fliplr(Err2Neg)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,Psth_Unrew,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

r1 = 0;
g1 = 0;
b1 = 255;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos, fliplr(Err1Neg)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,Psth_Rew,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
xlabel('Time (s)', 'FontSize', 10);
ylabel('Z-score', 'FontSize', 10);
line('XData', [(-1 * sec_preEvent),sec_postEvent], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [-1, 2], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
if c == 1
    legend ([h1 h2], 'Reward', 'No reward','orientation', 'vertical', 'Location', 'NorthWest');
    legend BOXOFF;
end
xlim ([(-1 * sec_preEvent),sec_postEvent]);
% ylim ([ymin, ymax]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 10, ...
  'LineWidth'   , 2         );
title(Photometry_zScore(c).subjectID);


%% Create PSTH to look at differential responses for early/late trials - photometry data

% Define number of trials to average across for early/late comparison
Num_early_vs_late_trials = 5;

% Rewarded entries
nRewLicks_ts = size(RewLicks_ts,1);
PsthArray_RL_early = PsthArray_Rew(1:Num_early_vs_late_trials, :);
PsthArray_RL_late = PsthArray_Rew(nRewLicks_ts-Num_early_vs_late_trials+1:end, :);

Psth_RL_early = nanmean(PsthArray_RL_early);
err_RL_early = (nanstd(PsthArray_RL_early))/sqrt(size(PsthArray_RL_early,1));

Psth_RL_late = nanmean(PsthArray_RL_late);
err_RL_late = (nanstd(PsthArray_RL_late))/sqrt(size(PsthArray_RL_late,1));

% Unrewarded entries
nUnrewLicks_ts = size(UnrewLicks_ts,1);
PsthArray_URL_early = PsthArray_Unrew(1:Num_early_vs_late_trials, :);
PsthArray_URL_late = PsthArray_Unrew(nUnrewLicks_ts-Num_early_vs_late_trials+1:end, :);

Psth_URL_early = nanmean(PsthArray_URL_early);
err_URL_early = (nanstd(PsthArray_URL_early))/sqrt(size(PsthArray_URL_early,1));

Psth_URL_late = nanmean(PsthArray_URL_late);
err_URL_late = (nanstd(PsthArray_URL_late))/sqrt(size(PsthArray_URL_late,1));


%% Plot early and late rewarded licks - photometry data
 
nTsPrev = round (sec_preEvent * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (sec_postEvent * samplingRate);
totalTs = nTsPrev + nTsPost;                     % sets time axis for upcoming graph based on these values
increment = (sec_preEvent + sec_postEvent) / totalTs;
timeAxis = (-1 * sec_preEvent) : increment : sec_postEvent;
ntimeAxis = size (timeAxis, 2);

Err1Pos = Psth_RL_early + err_RL_early;
Err1Neg = Psth_RL_early - err_RL_early;
Err2Pos = Psth_RL_late + err_RL_late;
Err2Neg = Psth_RL_late - err_RL_late;

if c== 1
   batch_RL_earlylate = figure;
end
figure (batch_RL_earlylate);
subplot (subject_plot_rows, subject_plot_cols, c);
hold on

% Change RGB values to adjust color, change 'o' value to adjust opacity for shaded error fills
o = 0.75;

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[Err2Pos(1, 1:ntimeAxis), fliplr(Err2Neg(1, 1:ntimeAxis))], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,Psth_RL_late(1, 1:ntimeAxis),'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3, 'LineStyle', ':');

o = 0.5;

r1 = 0;
g1 = 0;
b1 = 255;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos(1, 1:ntimeAxis), fliplr(Err1Neg(1, 1:ntimeAxis))], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,Psth_RL_early(1, 1:ntimeAxis),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
xlabel('Time (s)', 'FontSize', 10);
ylabel('Z-score', 'FontSize', 10);
line('XData', [(-1 * sec_preEvent),sec_postEvent], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [-1, 2], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
if c == 1
    legend ([h1 h2], 'Reward (early)', 'Reward (late)','orientation', 'vertical', 'Location', 'NorthWest');
    legend BOXOFF;
end
xlim ([(-1 * sec_preEvent),sec_postEvent]);
% ylim ([ymin, ymax]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 10, ...
  'LineWidth'   , 2         );

title(Photometry_zScore(c).subjectID);


%% Plot early and late unrewarded licks - photometry data

nTsPrev = round (sec_preEvent * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (sec_postEvent * samplingRate);
totalTs = nTsPrev + nTsPost;                     % sets time axis for upcoming graph based on these values
increment = (sec_preEvent + sec_postEvent) / totalTs;
timeAxis = (-1 * sec_preEvent) : increment : sec_postEvent;
ntimeAxis = size (timeAxis, 2);

Err1Pos = Psth_URL_early + err_URL_early;
Err1Neg = Psth_URL_early - err_URL_early;
Err2Pos = Psth_URL_late + err_URL_late;
Err2Neg = Psth_URL_late - err_URL_late;

if c== 1
   batch_URL_earlylate = figure;
end
figure (batch_URL_earlylate);
subplot (subject_plot_rows, subject_plot_cols, c);
hold on

% Change RGB values to adjust color, change 'o' value to adjust opacity for shaded error fills
o = 0.5;

r2 = 150;
g2 = 150;
b2 = 150;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[Err2Pos(1, 1:ntimeAxis), fliplr(Err2Neg(1, 1:ntimeAxis))], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,Psth_URL_late(1, 1:ntimeAxis),'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3, 'LineStyle', ':');

r1 = 75;
g1 = 75;
b1 = 75;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos(1, 1:ntimeAxis), fliplr(Err1Neg(1, 1:ntimeAxis))], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,Psth_URL_early(1, 1:ntimeAxis),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
xlabel('Time (s)', 'FontSize', 10);
ylabel('Z-score', 'FontSize', 10);
line('XData', [(-1 * sec_preEvent),sec_postEvent], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [-1, 2], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
if c == 1
    legend ([h1 h2], 'No reward (early)', 'No reward (late)','orientation', 'vertical', 'Location', 'NorthWest');
    legend BOXOFF;
end
xlim ([(-1 * sec_preEvent),sec_postEvent]);
% ylim ([ymin, ymax]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 10, ...
  'LineWidth'   , 2         );
title(Photometry_zScore(c).subjectID);


%% Store individual data in group structure arrays

% Behavioral data timestamps and PSTH
Behavior_ts(c).subjectID = Photometry_zScore(c).subjectID;
Behavior_ts(c).original_datafile = Photometry_zScore(c).original_datafile;

Behavior_ts(c).Licks_ts = Licks_ts;
Behavior_ts(c).Rewards_ts = Rewards_ts;
Behavior_ts(c).RewardedLicks_ts = RewLicks_ts;
Behavior_ts(c).UnrewardedLicks_ts = UnrewLicks_ts;

% Behavior lick rate PSTH + error (SEM) 
Behavior_PSTH(c).subjectID = Photometry_zScore(c).subjectID;
Behavior_PSTH(c).original_datafile = Photometry_zScore(c).original_datafile;

Behavior_PSTH(c).rewLick_PSTH = beh_rewLick_PSTH;
Behavior_PSTH(c).rewLick_err = beh_rewLick_err;
Behavior_PSTH(c).UnrewLick_PSTH = beh_UnrewLick_PSTH;
Behavior_PSTH(c).UnrewLick_err = beh_UnrewLick_err;

% Photometry PSTH and error
Photometry_PSTH(c).subjectID = Photometry_zScore(c).subjectID;
Photometry_PSTH(c).original_datafile = Photometry_zScore(c).original_datafile;
Photometry_PSTH(c).photometry_samplingrate = Photometry_zScore(c).photometry_samplingrate;

Photometry_PSTH(c).RewardLick_PSTH = Psth_Rew;
Photometry_PSTH(c).RewardLick_err = err_Rew;
Photometry_PSTH(c).UnrewardLick_PSTH = Psth_Unrew;
Photometry_PSTH(c).UnrewardLick_err = err_Unrew;

Photometry_PSTH(c).early_RewLick_PSTH = Psth_RL_early;
Photometry_PSTH(c).early_RewLick_err = err_RL_early;
Photometry_PSTH(c).late_RewLick_PSTH = Psth_RL_late;
Photometry_PSTH(c).late_RewLick_err = err_RL_late;

Photometry_PSTH(c).early_UnrewLick_PSTH = Psth_URL_early;
Photometry_PSTH(c).early_UnrewLick_err = err_URL_early;
Photometry_PSTH(c).late_UnrewLick_PSTH = Psth_URL_late;
Photometry_PSTH(c).late_UnrewLick_err = err_URL_late;

clearvars -except c Behavior* Photometry* batch* figname* subject_plot_rows subject_plot_cols behavior_files Num_early_vs_late_trials sec_preEvent sec_postEvent beh_time_resolution output_datafile_name;

end

%% Save figures and data using specified filenames

savefig (batch_licks_beh, figname1);
savefig (batch_licks, figname2);
savefig (batch_RL_earlylate, figname3);
savefig (batch_URL_earlylate, figname4);

clearvars -except Behavior* Photometry*  Num_early_vs_late_trials sec_preEvent sec_postEvent beh_time_resolution output_datafile_name;

save (output_datafile_name);
