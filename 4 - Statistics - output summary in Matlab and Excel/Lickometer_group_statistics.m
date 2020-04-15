%% Statistically analyze pre vs. post event data for behavior and neural activity in a lickometer task

% ** What you can expect this script to do:
%  Compare pre/post event behavior or neural activity using paired t-tests
%  Print summary results in excel spreadsheet

%  ** Key output data is stored in a structure array "Statistics" with 9 fields:
%       Name - Description of comparison
%       Comparitor1 - data1
%       Comparitor2 - data2
%       Hypothesis - 1 if null hypothesis (comparitors come from distributions with equal means) can be rejected at the 5% level, 0 if not
%       pValue - obvious
%       Mean1 - group mean for the first comparitor
%       Mean2 - group mean for the second comparitor
%       Values1 - comparitor 1 individual data used to generate group average
%       Values2 - comparitor 2 individual data used to generate group average

% ** What you need to customize:
%       Import/export data filenames
%       Analysis window

% ** Assumptions:
%   Raw data have been previously extracted from TDT data files and photometry signal has been fitted and z-score normalized using TDT_photometry_data_extraction_script
%   Data has been previously analyzed to generate subject-level PSTH for behavior and photometry with Lickometer_batch_analysis_github script


%% Prepare workspace

clear all
close all

%% Specify input and output filenames and time window for statistical analyses

input_datafile_name = 'Lickometer_group_analyzed_data';
stats_output_datafile_name = 'Lickometer_statistics';

analysis_window = 2;     % Pre/post event time window (seconds) for statistical comparison

%% Load previously analyzed individual data

x = 0;  % Initialize counter variable for statistics array

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
        
% Report an error if stats analysis time window exceeds original analysis
if analysis_window > originalSecPrev
    error ('The statistical analysis window exceeds the original analysis window for individual data.');
end

if analysis_window > originalSecPost 
    error ('The statistical analysis window exceeds the original analysis window for individual data.');
end


%% Statistics for behavior

% Extract individual PSTH vectors from structure array
tempCell = {Behavior_PSTH.rewLick_PSTH}';
group_behavior_rewLick_PSTH = cell2mat(tempCell);

tempCell = {Behavior_PSTH.UnrewLick_PSTH}';
group_behavior_UnrewLick_PSTH = cell2mat(tempCell);

% Determine bins to analyze
Beh_bins = analysis_window/beh_time_resolution;
Beh_zero = (originalSecPrev/beh_time_resolution)+1;
if Beh_bins == Beh_zero
    Beh_bins = Beh_bins-1;
end
Beh_totalBins = 1:size(group_behavior_rewLick_PSTH, 2);
Beh_preBins = Beh_totalBins(1, (Beh_zero-Beh_bins):1:(Beh_zero-1));
Beh_postBins = Beh_totalBins(1, Beh_zero:1:(Beh_zero+Beh_bins-1));

% Stats for rewarded licks
for i = 1:nSubjects
    Beh_Rew_Pre(i, 1) = mean(group_behavior_rewLick_PSTH(i, Beh_preBins));
    Beh_Rew_Post(i, 1) = mean(group_behavior_rewLick_PSTH(i, Beh_postBins));
end

Beh_Rew_ppDiff = Beh_Rew_Post - Beh_Rew_Pre;
Mean_Beh_Rew_ppDiff = mean(Beh_Rew_ppDiff);

Mean_Beh_Rew_Pre = mean(Beh_Rew_Pre);
Mean_Beh_Rew_Post = mean(Beh_Rew_Post);
[Beh_Rew_h, Beh_Rew_p] = ttest(Beh_Rew_Pre, Beh_Rew_Post);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Lick rate: rewarded lick';
Statistics(x).Comparitor1 = 'Lick rate: pre-rewarded lick';
Statistics(x).Comparitor2 = 'Lick rate: post-rewarded lick';
Statistics(x).Hypothesis = Beh_Rew_h;
Statistics(x).pValue = Beh_Rew_p;
Statistics(x).Mean1 = Mean_Beh_Rew_Pre;
Statistics(x).Mean2 = Mean_Beh_Rew_Post;
Statistics(x).Values1 = Beh_Rew_Pre;
Statistics(x).Values2 = Beh_Rew_Post;


% Stats for unrewarded licks
for i = 1:nSubjects
    Beh_NoRew_Pre(i, 1) = mean(group_behavior_UnrewLick_PSTH(i, Beh_preBins));
    Beh_NoRew_Post(i, 1) = mean(group_behavior_UnrewLick_PSTH(i, Beh_postBins));
end

Beh_NoRew_ppDiff = Beh_NoRew_Post - Beh_NoRew_Pre;
Mean_Beh_NoRew_ppDiff = mean(Beh_NoRew_ppDiff);

Mean_Beh_NoRew_Pre = mean(Beh_NoRew_Pre);
Mean_Beh_NoRew_Post = mean(Beh_NoRew_Post);
[Beh_NoRew_h, Beh_NoRew_p] = ttest(Beh_NoRew_Pre, Beh_NoRew_Post);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Lick rate: unrewarded lick';
Statistics(x).Comparitor1 = 'Lick rate: pre-unrewarded lick';
Statistics(x).Comparitor2 = 'Lick rate: post-unrewarded lick';
Statistics(x).Hypothesis = Beh_NoRew_h;
Statistics(x).pValue = Beh_NoRew_p;
Statistics(x).Mean1 = Mean_Beh_NoRew_Pre;
Statistics(x).Mean2 = Mean_Beh_NoRew_Post;
Statistics(x).Values1 = Beh_NoRew_Pre;
Statistics(x).Values2 = Beh_NoRew_Post;


% Compare conditions
[Beh_ppDiff_h, Beh_ppDiff_p] = ttest(Beh_NoRew_ppDiff, Beh_Rew_ppDiff);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Lick rate: rewarded vs. unrewarded';
Statistics(x).Comparitor1 = 'Lick rate: rewarded rate (post-pre)';
Statistics(x).Comparitor2 = 'Lick rate: unrewarded rate (post-pre)';
Statistics(x).Hypothesis = Beh_ppDiff_h;
Statistics(x).pValue = Beh_ppDiff_p;
Statistics(x).Mean1 = Mean_Beh_Rew_ppDiff;
Statistics(x).Mean2 = Mean_Beh_NoRew_ppDiff;
Statistics(x).Values1 = Beh_Rew_ppDiff;
Statistics(x).Values2 = Beh_NoRew_ppDiff;


%% Statistics for photometry (all trials reward/no reward)

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.RewardLick_PSTH}';
group_photometry_RewardLick_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.UnrewardLick_PSTH}';
group_photometry_UnrewardLick_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% Determine bins to analyze
interval = 1/samplingRate;
bins = round(analysis_window/interval);
zero = round(originalSecPrev/interval);
if bins == zero
    bins = bins-1;
end
totalBins = 1:size(group_photometry_RewardLick_PSTH, 2);
preBins = totalBins(1, (zero-bins):1:(zero-1));
postBins = totalBins(1, zero:1:(zero+bins-1));

% Stats for rewarded licks - photometry
for i = 1:nSubjects 
    Rew_Pre(i, 1) = mean(group_photometry_RewardLick_PSTH(i, preBins));
    Rew_Post(i, 1) = mean(group_photometry_RewardLick_PSTH(i, postBins));
end

Rew_ppDiff = Rew_Post - Rew_Pre;
Mean_Rew_ppDiff = mean(Rew_ppDiff);

Mean_Rew_Pre = mean(Rew_Pre);
Mean_Rew_Post = mean(Rew_Post);
[Rew_h, Rew_p] = ttest(Rew_Pre, Rew_Post);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: rewarded lick';
Statistics(x).Comparitor1 = 'Photometry: pre-rewarded lick';
Statistics(x).Comparitor2 = 'Photometry: post-rewarded lick';
Statistics(x).Hypothesis = Rew_h;
Statistics(x).pValue = Rew_p;
Statistics(x).Mean1 = Mean_Rew_Pre;
Statistics(x).Mean2 = Mean_Rew_Post;
Statistics(x).Values1 = Rew_Pre;
Statistics(x).Values2 = Rew_Post;


% Stats for unrewarded licks
for i = 1:nSubjects
    NoRew_Pre(i, 1) = mean(group_photometry_UnrewardLick_PSTH(i, preBins));
    NoRew_Post(i, 1) = mean(group_photometry_UnrewardLick_PSTH(i, postBins));
end

NoRew_ppDiff = NoRew_Post - NoRew_Pre;
Mean_NoRew_ppDiff = mean(NoRew_ppDiff);

Mean_NoRew_Pre = mean(NoRew_Pre);
Mean_NoRew_Post = mean(NoRew_Post);
[NoRew_h, NoRew_p] = ttest(NoRew_Pre, NoRew_Post);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: unrewarded lick';
Statistics(x).Comparitor1 = 'Photometry: pre-unrewarded lick';
Statistics(x).Comparitor2 = 'Photometry: post-unrewarded lick';
Statistics(x).Hypothesis = NoRew_h;
Statistics(x).pValue = NoRew_p;
Statistics(x).Mean1 = Mean_NoRew_Pre;
Statistics(x).Mean2 = Mean_NoRew_Post;
Statistics(x).Values1 = NoRew_Pre;
Statistics(x).Values2 = NoRew_Post;

% Compare conditions
[ppDiff_h, ppDiff_p] = ttest(NoRew_ppDiff, Rew_ppDiff);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: rewarded vs. unrewarded lick';
Statistics(x).Comparitor1 = 'Photometry: rewarded lick (post-pre)';
Statistics(x).Comparitor2 = 'Photometry: unrewarded lick (post-pre)';
Statistics(x).Hypothesis = ppDiff_h;
Statistics(x).pValue = ppDiff_p;
Statistics(x).Mean1 = Mean_Rew_ppDiff;
Statistics(x).Mean2 = Mean_NoRew_ppDiff;
Statistics(x).Values1 = Rew_ppDiff;
Statistics(x).Values2 = NoRew_ppDiff;

%% Statistics for photometry (Early vs. late reward) - calculate average value across all intervals in defined time bin pre/post event

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.early_RewLick_PSTH}';
group_photometry_early_RewLick_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.late_RewLick_PSTH}';
group_photometry_late_RewLick_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.early_UnrewLick_PSTH}';
group_photometry_early_UnrewLick_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.late_UnrewLick_PSTH}';
group_photometry_late_UnrewLick_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% Determine bins to analyze
interval = 1/samplingRate;
bins = round(analysis_window/interval);
zero = round(originalSecPrev/interval);
if bins == zero
    bins = bins-1;
end
totalBins = 1:size(group_photometry_early_RewLick_PSTH, 2);
preBins = totalBins(1, (zero-bins):1:(zero-1));
postBins = totalBins(1, zero:1:(zero+bins-1));

% Stats for early and late rewarded licks - photometry
for i = 1:nSubjects
    Rew_EarlyPre(i, 1) = mean(group_photometry_early_RewLick_PSTH(i, preBins));
    Rew_EarlyPost(i, 1) = mean(group_photometry_early_RewLick_PSTH(i, postBins));
    Rew_LatePre(i, 1) = mean(group_photometry_late_RewLick_PSTH(i, preBins));
    Rew_LatePost(i, 1) = mean(group_photometry_late_RewLick_PSTH(i, postBins));
end

Rew_Early_ppDiff = Rew_EarlyPost - Rew_EarlyPre;
Mean_Rew_Early_ppDiff = mean(Rew_Early_ppDiff);

Rew_Late_ppDiff = Rew_LatePost - Rew_LatePre;
Mean_Rew_Late_ppDiff = mean(Rew_Late_ppDiff);

% Compare conditions
[Rew_EL_ppDiff_h, Rew_EL_ppDiff_p] = ttest(Rew_Early_ppDiff, Rew_Late_ppDiff);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: early vs. late rewarded trials';
Statistics(x).Comparitor1 = 'early rewarded trials (pre-post difference)';
Statistics(x).Comparitor2 = 'late rewarded trials (pre-post difference)';
Statistics(x).Hypothesis = Rew_EL_ppDiff_h;
Statistics(x).pValue = Rew_EL_ppDiff_p;
Statistics(x).Mean1 = Mean_Rew_Early_ppDiff;
Statistics(x).Mean2 = Mean_Rew_Late_ppDiff;
Statistics(x).Values1 = Rew_Early_ppDiff;
Statistics(x).Values2 = Rew_Late_ppDiff;

%Stats for early and late unrewarded licks - photometry
for i = 1:nSubjects
    NoRew_EarlyPre(i, 1) = mean(group_photometry_early_UnrewLick_PSTH(i, preBins));
    NoRew_EarlyPost(i, 1) = mean(group_photometry_early_UnrewLick_PSTH(i, postBins));
    NoRew_LatePre(i, 1) = mean(group_photometry_late_UnrewLick_PSTH(i, preBins));
    NoRew_LatePost(i, 1) = mean(group_photometry_late_UnrewLick_PSTH(i, postBins));
end

NoRew_Early_ppDiff = NoRew_EarlyPost - NoRew_EarlyPre;
Mean_NoRew_Early_ppDiff = mean(NoRew_Early_ppDiff);

NoRew_Late_ppDiff = NoRew_LatePost - NoRew_LatePre;
Mean_NoRew_Late_ppDiff = mean(NoRew_Late_ppDiff);

% Compare conditions
[NoRew_EL_ppDiff_h, NoRew_EL_ppDiff_p] = ttest(NoRew_Early_ppDiff, NoRew_Late_ppDiff);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: early vs. late unrewarded trials';
Statistics(x).Comparitor1 = 'early unrewarded trials (pre-post difference)';
Statistics(x).Comparitor2 = 'late unrewarded trials (pre-post difference)';
Statistics(x).Hypothesis = NoRew_EL_ppDiff_h;
Statistics(x).pValue = NoRew_EL_ppDiff_p;
Statistics(x).Mean1 = Mean_NoRew_Early_ppDiff;
Statistics(x).Mean2 = Mean_NoRew_Late_ppDiff;
Statistics(x).Values1 = NoRew_Early_ppDiff;
Statistics(x).Values2 = NoRew_Late_ppDiff;


%% Save data

% Print summary to excel
Excelprint = Statistics;
Excelprint = rmfield(Excelprint, 'Values1');
Excelprint = rmfield(Excelprint, 'Values2');
writetable(struct2table(Excelprint), 'Lickometer Stats Summary.xlsx')

clearvars -except Statistics Num_early_vs_late_trials analysis_window Behavior_PSTH Photometry_PSTH input_datafile_name stats_output_datafile_name
save (stats_output_datafile_name);
