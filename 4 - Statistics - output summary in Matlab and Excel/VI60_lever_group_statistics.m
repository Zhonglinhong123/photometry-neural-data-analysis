%% Statistically analyze pre vs. post event data for neural activity in VI60 lever press task

% ** What you can expect this script to do:
%  Compare pre/post event behavior or neural activity using paired t-tests
%  Print summary results in excel spreadsheet

%  ** Key output data is stored in a structure array "Statistics" with 9 fields:
%       Name - Description of comparison
%       Comparitor1 - name of data1
%       Comparitor2 - name of data2
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
%   Data has been previously analyzed to generate subject-level PSTH for behavior and photometry with VI60_lever_batch_analysis script


%% Prepare workspace

clear all
close all

%% Specify input and output filenames and time window for statistical analyses

input_datafile_name = 'VI60_group_analyzed_data';
stats_output_datafile_name = 'VI60_lever_statistics';

analysis_window = 5;     % Pre/post event time window (seconds) for statistical comparison

%% Load previously analyzed individual data

x = 0;  % Initialize counter variable for statistics array

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

% Report an error if stats analysis time window exceeds original analysis
if analysis_window > originalSecPrev
    error ('The statistical analysis window exceeds the original analysis window for individual data.');
end

if analysis_window > originalSecPost 
    error ('The statistical analysis window exceeds the original analysis window for individual data.');
end

%% Statistics for photometry (rewarded vs unrewarded port entries)

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.RewardEntry_PSTH}';
group_photometry_RewardEntry_PSTH = cell2mat(tempCell);

tempCell = {Photometry_PSTH.UnrewardEntry_PSTH}';
group_photometry_UnrewardEntry_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% Determine bins to analyze
interval = 1/samplingRate;
bins = round(analysis_window/interval);
zero = round(originalSecPrev/interval);
if bins == zero
    bins = bins-1;
end
totalBins = 1:size(group_photometry_RewardEntry_PSTH, 2);
preBins = totalBins(1, (zero-bins):1:(zero-1));
postBins = totalBins(1, zero:1:(zero+bins-1));

% Stats for rewarded entries - photometry
for i = 1:nSubjects 
    Rew_Pre(i, 1) = mean(group_photometry_RewardEntry_PSTH(i, preBins));
    Rew_Post(i, 1) = mean(group_photometry_RewardEntry_PSTH(i, postBins));
end

Rew_ppDiff = Rew_Post - Rew_Pre;
Mean_Rew_ppDiff = mean(Rew_ppDiff);

Mean_Rew_Pre = mean(Rew_Pre);
Mean_Rew_Post = mean(Rew_Post);
[Rew_h, Rew_p] = ttest(Rew_Pre, Rew_Post);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: rewarded port entries (pre vs. post)';
Statistics(x).Comparitor1 = 'Pre-rewarded entry';
Statistics(x).Comparitor2 = 'Post-rewarded entry';
Statistics(x).Hypothesis = Rew_h;
Statistics(x).pValue = Rew_p;
Statistics(x).Mean1 = Mean_Rew_Pre;
Statistics(x).Mean2 = Mean_Rew_Post;
Statistics(x).Values1 = Rew_Pre;
Statistics(x).Values2 = Rew_Post;


% Stats for unrewarded entries
for i = 1:nSubjects
    NoRew_Pre(i, 1) = mean(group_photometry_UnrewardEntry_PSTH(i, preBins));
    NoRew_Post(i, 1) = mean(group_photometry_UnrewardEntry_PSTH(i, postBins));
end

NoRew_ppDiff = NoRew_Post - NoRew_Pre;
Mean_NoRew_ppDiff = mean(NoRew_ppDiff);

Mean_NoRew_Pre = mean(NoRew_Pre);
Mean_NoRew_Post = mean(NoRew_Post);
[NoRew_h, NoRew_p] = ttest(NoRew_Pre, NoRew_Post);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: unrewarded port entries (pre vs. post)';
Statistics(x).Comparitor1 = 'Pre-unrewarded entry';
Statistics(x).Comparitor2 = 'Post-unrewarded entry';
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
Statistics(x).Name = 'Photometry: rewarded vs. unrewarded port entry';
Statistics(x).Comparitor1 = 'Rewarded entry (post-pre)';
Statistics(x).Comparitor2 = 'Unrewarded entry (post-pre)';
Statistics(x).Hypothesis = ppDiff_h;
Statistics(x).pValue = ppDiff_p;
Statistics(x).Mean1 = Mean_Rew_ppDiff;
Statistics(x).Mean2 = Mean_NoRew_ppDiff;
Statistics(x).Values1 = Rew_ppDiff;
Statistics(x).Values2 = NoRew_ppDiff;

%% Statistics for photometry (active lever press) - calculate average value across all intervals in defined time bin pre/post event

% Extract individual PSTH vectors from structure array
tempCell = {Photometry_PSTH.Active_Lever_PSTH}';
group_photometry_Active_Lever_PSTH = cell2mat(tempCell);

samplingRate = Photometry_PSTH.photometry_samplingrate;

% Determine bins to analyze
interval = 1/samplingRate;
bins = round(analysis_window/interval);
zero = round(originalSecPrev/interval);
if bins == zero
    bins = bins-1;
end
totalBins = 1:size(group_photometry_Active_Lever_PSTH, 2);
preBins = totalBins(1, (zero-bins):1:(zero-1));
postBins = totalBins(1, zero:1:(zero+bins-1));

% Stats for lever press pre-post
for i = 1:nSubjects
    Lever_Pre(i, 1) = mean(group_photometry_Active_Lever_PSTH(i, preBins));
    Lever_Post(i, 1) = mean(group_photometry_Active_Lever_PSTH(i, postBins));
end

Mean_Lever_Pre = mean(Lever_Pre);
Mean_Lever_Post = mean(Lever_Post);
[Lever_h, Lever_p] = ttest(Lever_Pre, Lever_Post);

% Record results in structure array
x = x+1;
Statistics(x).Name = 'Photometry: active lever (pre vs. post)';
Statistics(x).Comparitor1 = 'Pre lever press)';
Statistics(x).Comparitor2 = 'Post lever press';
Statistics(x).Hypothesis = Lever_h;
Statistics(x).pValue = Lever_p;
Statistics(x).Mean1 = Mean_Lever_Pre;
Statistics(x).Mean2 = Mean_Lever_Post;
Statistics(x).Values1 = Lever_Pre;
Statistics(x).Values2 = Lever_Post;


%% Save data

% Print summary to excel
Excelprint = Statistics;
Excelprint = rmfield(Excelprint, 'Values1');
Excelprint = rmfield(Excelprint, 'Values2');
writetable(struct2table(Excelprint), 'VI60 Stats Summary.xlsx')

clearvars -except Statistics analysis_window Photometry_PSTH input_datafile_name stats_output_datafile_name
save (stats_output_datafile_name);
