%% TDT photometry data extraction and pre-processing %%%%%%%%%%%%%%%%%%%%%%

%  ** What you can expect this script to do:
%   Control signal is fitted to calcium dependent signal and subtracted to correct for motion artifacts and bleaching
%   Fitted signal is z-score normalized to facilitate comparisons across subjects
%   Graphs are generated for each subject of raw, fitted and z-scored data for QC check
%   Key output data stored in "Photometry_zScore" structure arraycontaining the following fields:
%       (1) Subject IDs, (2) Raw datafile names, (3) Sampling rate, 
%       (4) Fitted/normalized photometry data, (5) Photometry timestamps, (6+) TTL timestamps for external events recorded via TDT system

%  ** What you need to customize:
%       Import/export data filenames
%       TDT data sources(photometry data + TTL inputs to sync behavior/stimuli)

%  ** Assumptions:
%   2-color input: 470 nm (calcium dependent signal) and 405 nm (calcium independent control)
%   Requires additional standalone functions: tdt2mat, controlFit, deltaFF (credit: Tom Davidson)

%% Prepare workspace

clear all
close all

%% Specify subject identifiers, import/export file names, TDT tank name, TDT data sources

% Text identifier for each mouse
subject_IDs = {'FP1' ;'FP4' ;'FP11';'FP12'};                                                                       
nSubjects = numel(subject_IDs);

% Specify raw data import file names
indiv_raw_datafiles = {'FP1_VI60_100716'; 'FP4_VI60_100716'; 'FP11_VI60_100716'; 'FP12_VI60_100716a'};   

% Specify export filenames for individual signal QC graphs
indiv_signal_QC_graphs = {'FP1_VI60_100716_signal'; 'FP4_VI60_100716_signal'; 'FP11_VI60_100716_signal'; 'FP12_VI60_100716_signal'};

% Specify export file name for z-score normalized group photometry 
export_datafile_name = 'VI60_group_raw_data.mat';

% Point to folder where data is
data_import_folder = '/Volumes/Fake Server 2/Fiber photometry/Photometry data (all)/FP 1-18 Sept-Nov 2016'; 

% Specify TDT data tank name
TDT_datatank_name = 'Liz';

% Specify TDT data sources
% LMag = photometry signal; Din0 = optional TTL input for syncing behavioral events or stimuli (add more if needed)
TDT_data_sources = {'LMag' 'Din0' 'Din1' 'Din2' 'Din4' };   % LMag must be listed first. Specify additional TTL inputs (digital inputs, Din0...DinX) as needed - 3 places to modify code below if more inputs are added                                  
                                                                                                    
% Report an error in case of mismatch between subject IDs and filenames to ensure accurate labeling
if numel(indiv_raw_datafiles) ~= nSubjects | numel(indiv_signal_QC_graphs) ~= nSubjects
    error ('The number of subjects does not match the number of associated filenames.')
end

%%  Initiate analysis loop for each subject

for c = 1:nSubjects

%% Extract photometry data from TDT data tank

% Set up variables to extract
blockname = char(indiv_raw_datafiles(c));
figname1 = char(indiv_signal_QC_graphs(c));  %figname for three signal graphs

% extract
for k = 1:numel(TDT_data_sources)
  storename = TDT_data_sources {k};
  S{k} = tdt2mat(data_import_folder, TDT_datatank_name, blockname, storename);
end

%% Massage data and get time stamps

% LMag = photometry data
LMag = S{1};

% Din0...Din4 = TTL inputs (optional) 
Din0 = S{2};        % Import data from additional TTL inputs as needed - modify in 2 other places below
Din1 = S{3};
Din2 = S{4};  
Din4 = S{5};

% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 469nm, channel 2 = 405nm
chani1 = LMag.channels==1;
chani2 = LMag.channels==2;

% Get LMag data as a vector (repeat for each channel)
dat1 = LMag.data(chani1,:);
dat1 = reshape(dat1', [],1); % unwrap data from m x 256 array
dat2 = LMag.data(chani2,:);
dat2 = reshape(dat2', [],1); % unwrap data from m x 256 array

% Get LMag sampling rate
samplingRate = LMag.sampling_rate;

% Get LMag timestamps (use chani1 - timestamps should be the same for all LMag channels
ts = LMag.timestamps(chani1);
t_rec_start = ts(1);
ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
ts = reshape(ts',[],1);

% Get behavior (TTL input) timestamps
TTL1_timestamps = Din0.timestamps - t_rec_start;     % Include data from additional TTL inputs as needed
TTL2_timestamps = Din1.timestamps - t_rec_start;
TTL3_timestamps = Din2.timestamps - t_rec_start;
%TTL4_timestamps = Din3.timestamps - t_rec_start;
TTL5_timestamps = Din4.timestamps - t_rec_start;

%% Visualize photometry data (calcium-dependent signal and control) for QC check

x = figure;
subplot (4, 1, 1);
hold on; 
h1 = plot(ts, dat1(:), 'b');
h2 = plot(ts, dat2(:), 'r');
legend ([h1 h2], 'Raw signal', 'Raw control', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
xlabel('Time (s)', 'FontSize', 10);
ylabel('Raw data', 'FontSize', 10);
title(subject_IDs(c));

% Smooth dat1 and dat2, fit dat2 (control) to dat1 (signal)

dat1 = filtfilt(ones(1,100)/100,1, dat1);
dat2 = filtfilt(ones(1,100)/100,1, dat2);

[controlFit] = controlFit (dat1, dat2);

subplot (4, 1, 2);
hold on;
h1 = plot (ts, dat1, 'b');
h2 = plot (ts, controlFit, 'r');
legend ([h1 h2], 'Raw signal', 'Fitted control', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
xlabel('Time (s)', 'FontSize', 10);
ylabel('Fitted data', 'FontSize', 10);

% Get delta F/F using controlFit function

[normDat] = deltaFF (dat1, controlFit);

subplot (4, 1, 3);
plot (ts, normDat, 'k');
xlabel('Time (s)', 'FontSize', 10);
ylabel('\DeltaF/F (%)', 'FontSize', 10);

% Calculate Zscore
meanDat = mean(normDat);
stdevDat = std(normDat);
zDat = (normDat-meanDat)/stdevDat;

subplot (4, 1, 4);
plot (ts, zDat, 'k');
xlabel('Time (s)', 'FontSize', 10);
ylabel('Zscore', 'FontSize', 10);

savefig (x, figname1);

% Store individual data to group structure array
Photometry_zScore(c).subjectID = subject_IDs(c);
Photometry_zScore(c).original_datafile = indiv_raw_datafiles(c);
Photometry_zScore(c).photometry_samplingrate = samplingRate;
Photometry_zScore(c).photometry_zScore_data = zDat;
Photometry_zScore(c).photometry_timestamps = ts;
Photometry_zScore(c).TTL_event1 = TTL1_timestamps;
Photometry_zScore(c).TTL_event2 = TTL2_timestamps;     % Store data from additional TTL inputs as needed
Photometry_zScore(c).TTL_event3 = TTL3_timestamps;
%Photometry_zScore(c).TTL_event4 = TTL4_timestamps;
Photometry_zScore(c).TTL_event5 = TTL5_timestamps;

clearvars -except Photometry_zScore c export_datafile_name subject_IDs indiv_raw_datafiles nSubjects indiv_signal_QC_graphs data_import_folder TDT_datatank_name TDT_data_sources

end

%% Save pre-process photometry data to group file for later analysis

save (export_datafile_name);