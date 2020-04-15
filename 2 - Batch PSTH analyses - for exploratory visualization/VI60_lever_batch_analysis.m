%% Create subject-level PSTH for neural activity in a VI60 lever press task for reward

% ** What you can expect this script to do:
%  Select timestamps for rewarded and unrewarded port entries and lever presses that are separated in time to avoid re-sampling the same data twice
%  Generate batch graphs of photometry responses as PSTH for QC and exploratory visualization of individual variability

%  Key output data stored in two structure arrays: 
 
%       Photometry_PSTH - Structure array with individual-level data in 9 fields:
%       (1) subjectID (2) original datafile (3) sampling rate
%       (4) Rewarded Entry PSTH and (5) SEM
%       (6) Unrewarded Entry PSTH and (7) SEM
%       (8) Active lever PSTH and (9) SEM

%       Behavior_ts - Timestamps for behavioral data. Structure array with individual-level data 9 fields:
%       (1) subjectID (2) original datafile
%       (3) All active lever timestamps (4) All port entry timestamps
%       (5) Reward timestamps (6) All port exit timestamps 
%       (7) Rewarded port entry timestamps (8) Unrewarded entry timestamps (processed subset)
%       (9) Active lever timestamps for graph (processed subset with no overlap)

% ** What you need to customize:
%       Import/export data filenames
%       Source of behavioral data - separate matlab file or from TTL inputs through TDT system (will default to TTL inputs through TDT system if behavior filenames are not supplied)
%       Specify which TTL inputs correspond to behavioral events

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

input_datafile_name = 'VI60_group_raw_data';
output_datafile_name = 'VI60_group_analyzed_data';
figname1 = 'VI60_Rew_v_NoRew_PortEntry_individual_data';
figname2 = 'VI60_active_lever_individual_data';
subject_plot_rows = 2;                                                      
subject_plot_cols = 2; 

% Optional - specify separate files with behavioral data - specify variables to import below
behavior_files = {};

% Specify PSTH time window
sec_preEvent = 10;
sec_postEvent = 10;

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
    % Include your own code to import behavioral data if this applies
else
% Option 2 - use this code if behavioral events are captured in the same TDT datafile as photometry signal
    % Get behavior (TTL input) timestamps
    ActLever_ts = Photometry_zScore(c).TTL_event1;
    PortEntry_ts = Photometry_zScore(c).TTL_event2;
    Reward_ts = Photometry_zScore(c).TTL_event3;
    PortExit_ts = Photometry_zScore(c).TTL_event5;
end

%% Load photometry data

zDat = Photometry_zScore(c).photometry_zScore_data;
samplingRate = Photometry_zScore(c).photometry_samplingrate;

%% PORT ENTRIES - quality control
% (1) remove entries without exits (extremely short, correct for hardware limitation); 
% (2) calculate length of all port entries, 
% (3) remove entries lasting less than (criterion) sec 

% Determine whether there are the same number of entries and exits; remove entries that have no corresponding exits (presumably so fast that exit not detected)
nPortEntry = size(PortEntry_ts, 1);
nPortExit = size(PortExit_ts, 1);
PortMismatch = nPortEntry - nPortExit;
if PortMismatch ~= 0;
    if PortEntry_ts(nPortEntry, 1) > PortExit_ts(nPortExit, 1);   %This detects if mouse was in port when session ended (i.e., last entry happened after last exit. If so, delete last entry to prevent code errors).
       PortEntry_ts (nPortEntry, 1) = NaN;
       PortEntry_ts = PortEntry_ts(all(~isnan(PortEntry_ts), 2),:);
       nPortEntry = size(PortEntry_ts, 1);
       PortMismatch = nPortEntry - nPortExit;
    end
    for iPort = 1:(nPortEntry-PortMismatch-1);
        thisEntry = PortEntry_ts(iPort);
        nextEntry = PortEntry_ts(iPort+1);
        thisExit = PortExit_ts(iPort);
        if (nextEntry - thisExit) <= 0;
           PortEntry_ts(iPort) = NaN;       %delete port entry that has no corresponding exit here
           PortEntry_ts = PortEntry_ts(all(~isnan(PortEntry_ts), 2),:);
        end
    end
end

% Calculate duration for each port entry/exit combination
nPortDuration = size(PortEntry_ts, 1);
PortDuration = NaN(nPortDuration, 1);
for iDur = 1:nPortDuration;
    thisEntry = PortEntry_ts(iDur);
    thisExit = PortExit_ts(iDur);
    PortDuration(iDur, 1) = (thisExit - thisEntry);
end

% Exclude port entries less than criterion sec
PortCriterion = 0.02;  
for iDur = 1:nPortDuration;   
    if PortDuration (iDur, 1) <= PortCriterion;
        PortDuration (iDur, 1) = NaN;
    end
end
PortLogical = ~isnan(PortDuration);
ShortPort = sum(isnan(PortDuration));
PortDuration = PortDuration(PortLogical);
PortEntry_ts = PortEntry_ts(PortLogical);
PortExit_ts = PortExit_ts(PortLogical);


%% PORT ENTRIES - separate rewarded and unrewarded port entries 

% Find first port entry following reward 
nReward_ts = size(Reward_ts, 1);
nPortEntry = size(PortEntry_ts, 1);
RewEntry_ts = NaN(nReward_ts, 1);
InPort = 0;
for iRew = 1:nReward_ts;
    thisRew = Reward_ts(iRew);
    for iPort = 1:nPortEntry;
        if thisRew >= PortEntry_ts(iPort, 1) & thisRew <= PortExit_ts(iPort, 1);
           RewEntry_ts(iRew,1) = thisRew;                       % rewarded entry = reward time if mouse is already in port
           InPort = InPort+1;
        end
    end
    if isnan(RewEntry_ts(iRew,1));
       thisEntry = find(PortEntry_ts>=thisRew, 1, 'first');
       if isempty (thisEntry) == 1;
           RewEntry_ts(iRew,1) = NaN;                           % delete reward if no entries occurred during or after
       else
           RewEntry_ts(iRew,1) = PortEntry_ts(thisEntry,1);     % rewarded entry = 1st PE after reward if mouse not in port when reward delivered
       end
    end
end
RewEntry_ts = RewEntry_ts(all(~isnan(RewEntry_ts), 2),:);

% Find unrewarded port entries (multistep calculation)

PreRewCriterion = 20;                       %must fall outside of criterion times before and after rewarded port entries
PostRewCriterion = 20;
UnrewEntry_ts = PortEntry_ts;
nUnrewEntry = size(UnrewEntry_ts, 1);
nRewEntry = size(RewEntry_ts, 1);
for iRew = 1:nRewEntry;                     %remove reward PE and surrounding entries
    thisRewEntry = RewEntry_ts(iRew, 1);
    for iUnRew = 1:nUnrewEntry;
        if UnrewEntry_ts(iUnRew, 1) > (thisRewEntry - PreRewCriterion) & UnrewEntry_ts(iUnRew, 1) < (thisRewEntry + PostRewCriterion);
           UnrewEntry_ts(iUnRew, 1) = NaN;
        end
    end
end
UnrewEntry_ts = UnrewEntry_ts(all(~isnan(UnrewEntry_ts), 2),:);     %Unrewarded PE times - BUT these may occur in bursts which are problematic for PSTH averaging



PostUnrewCriterion = 20;                    %remove unrewarded PE that occur within criterion sec of prev unrewarded PE (to avoid duplicating same neural data with slight time shift in PSTH)
nUnrewEntry = size(UnrewEntry_ts, 1);
for iUnRew = 1:nUnrewEntry;
    remainder = nUnrewEntry - iUnRew;
    for i = 1:remainder;
        if UnrewEntry_ts(iUnRew+i, 1) <= (UnrewEntry_ts(iUnRew, 1) + PostUnrewCriterion);
           UnrewEntry_ts(iUnRew+i, 1) = NaN;
        end
    end
end
UnrewEntry_ts = UnrewEntry_ts(all(~isnan(UnrewEntry_ts), 2),:); % This is heavily processed set of PE times that meet all criteria for being unrewarded and are separated in time by at least 15s


%% LEVER PRESSES - Select active lever presses that don't overlap during analysis window

lever_crit = 1; %set criterion for lever analysis window - exclude behavioral events happening criterion sec before/after
nLever_ts = size(ActLever_ts, 1);    
ActLever_no_overlap_ts = ActLever_ts;

for i = 1:nLever_ts;
    remainder = nLever_ts - i;
    for j = 1:remainder;
        if ActLever_no_overlap_ts(i+j, 1) <= (ActLever_no_overlap_ts(i, 1) + lever_crit);
           ActLever_no_overlap_ts(i+j, 1) = NaN;
        end
    end
end
ActLever_no_overlap_ts = ActLever_no_overlap_ts(all(~isnan(ActLever_no_overlap_ts), 2),:);


%% Make PSTH arrays of photometry data for rewarded and unrewarded licks, then average rows into a mean vector for plotting

% Convert seconds to TDT timestamps
nTsPrev = round (sec_preEvent * samplingRate);
nTsPost = round (sec_postEvent * samplingRate);


% Make PSTH for rewarded entries
nRewEntry = size(RewEntry_ts,1);
PsthArray_Rew = NaN(nRewEntry,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nRewEntry
    thisTime = RewEntry_ts(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_Rew(i,:) = processPhotDataRow_normDat(zDat, thisIndex, nTsPrev, nTsPost);
end
err_Rew = (nanstd(PsthArray_Rew))/sqrt(size(PsthArray_Rew,1));
Psth_Rew = nanmean(PsthArray_Rew);


% Make PSTH for unrewarded entries
nUnrewEntry = size(UnrewEntry_ts,1);
PsthArray_Unrew = NaN(nUnrewEntry,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nUnrewEntry
    thisTime = UnrewEntry_ts(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_Unrew(i,:) = processPhotDataRow_normDat(zDat, thisIndex, nTsPrev, nTsPost);
end
err_Unrew = (nanstd(PsthArray_Unrew))/sqrt(size(PsthArray_Unrew,1));
Psth_Unrew = nanmean(PsthArray_Unrew);


% Make PSTH for active lever presses
nActLever = size(ActLever_no_overlap_ts,1);
PsthArray_ActLever = NaN(nActLever,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nActLever
    thisTime = ActLever_no_overlap_ts(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_ActLever(i,:) = processPhotDataRow_normDat(zDat, thisIndex, nTsPrev, nTsPost);
end
err_Lever = (nanstd(PsthArray_ActLever))/sqrt(size(PsthArray_ActLever,1));
Psth_Lever = nanmean(PsthArray_ActLever);


%% Plot rewarded and unrewarded entries - photometry data

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
   batch_port_entries = figure;
end
figure (batch_port_entries);
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


%% Plot active lever responses - photometry data
 
nTsPrev = round (sec_preEvent * samplingRate);       % convert seconds to TDT timestamps
nTsPost = round (sec_postEvent * samplingRate);
totalTs = nTsPrev + nTsPost;                     % sets time axis for upcoming graph based on these values
increment = (sec_preEvent + sec_postEvent) / totalTs;
timeAxis = (-1 * sec_preEvent) : increment : sec_postEvent;
ntimeAxis = size (timeAxis, 2);

Err1Pos = Psth_Lever + err_Lever;
Err1Neg = Psth_Lever - err_Lever;


if c== 1
   batch_activeLever = figure;
end
figure (batch_activeLever);
subplot (subject_plot_rows, subject_plot_cols, c);
hold on

% Change RGB values to adjust color, change 'o' value to adjust opacity for shaded error fills

o = 0.5;

r1 = 24;
g1 = 86;
b1 = 9;
rgb1_o = opacity (o, r1, g1, b1);

fill([timeAxis, fliplr(timeAxis)],[Err1Pos(1, 1:ntimeAxis), fliplr(Err1Neg(1, 1:ntimeAxis))], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,Psth_Lever(1, 1:ntimeAxis),'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
xlabel('Time (s)', 'FontSize', 10);
ylabel('Z-score', 'FontSize', 10);
line('XData', [(-1 * sec_preEvent),sec_postEvent], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [-1, 2], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
if c == 1
    legend ([h1], 'Active lever','orientation', 'vertical', 'Location', 'NorthWest');
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

% Processed behavioral data timestamps
Behavior_ts(c).subjectID = Photometry_zScore(c).subjectID;
Behavior_ts(c).original_datafile = Photometry_zScore(c).original_datafile;

Behavior_ts(c).Active_Lever_ALL_ts = ActLever_ts;
Behavior_ts(c).PortEntry_ts = ActLever_ts;
Behavior_ts(c).Reward_ts = ActLever_ts;
Behavior_ts(c).PortExit_ts = ActLever_ts;
Behavior_ts(c).RewardedEntry_ts = RewEntry_ts;
Behavior_ts(c).UnrewardedEntry_ts = UnrewEntry_ts;
Behavior_ts(c).Active_Lever_forGraph_ts = ActLever_no_overlap_ts;

% Photometry PSTH and error
Photometry_PSTH(c).subjectID = Photometry_zScore(c).subjectID;
Photometry_PSTH(c).original_datafile = Photometry_zScore(c).original_datafile;
Photometry_PSTH(c).photometry_samplingrate = Photometry_zScore(c).photometry_samplingrate;

Photometry_PSTH(c).RewardEntry_PSTH = Psth_Rew;
Photometry_PSTH(c).RewardEntry_err = err_Rew;
Photometry_PSTH(c).UnrewardEntry_PSTH = Psth_Unrew;
Photometry_PSTH(c).UnrewardEntry_err = err_Unrew;

Photometry_PSTH(c).Active_Lever_PSTH = Psth_Lever;
Photometry_PSTH(c).Active_Lever_err = err_Lever;

clearvars -except c Behavior* Photometry* batch* figname* subject_plot_rows subject_plot_cols behavior_files sec_preEvent sec_postEvent beh_time_resolution output_datafile_name;

end

%% Save figures and data using specified filenames

savefig (batch_port_entries, figname1);
savefig (batch_activeLever, figname2);

clearvars -except Behavior* Photometry*  Num_early_vs_late_trials sec_preEvent sec_postEvent beh_time_resolution output_datafile_name;

save (output_datafile_name);
