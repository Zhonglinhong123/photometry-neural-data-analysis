# photometry-neural-data-analysis

Matlab scripts for batch import and analysis of photometry data from TDT-based systems.

Data processing is divided into four phases; behavior-specific files differ in behavioral data analysis specifics but share common strategy for processing photometry data and generating PSTHs.

Some statistical and behavioral analyses reported in Steinberg et al. Neuron 2020 are not included in this code as they were conducted in other programs.

***********************************************
1) Data extraction and preparation for analysis, QC check

What you can expect this script to do:
        Control signal is fitted to calcium dependent signal and subtracted to correct for motion artifacts and bleaching
        Fitted signal is z-score normalized to facilitate comparisons across subjects
        Graphs are generated for each subject of raw, fitted and z-scored data for QC check

Key output data stored in "Photometry_zScore" structure array containing the following fields:
       (1) Subject IDs, (2) Raw datafile names, (3) Sampling rate, 
       (4) Fitted/normalized photometry data, (5) Photometry timestamps, (6+) TTL timestamps for external events
       
What you need to customize:
       Input/output data filenames
       Specify TDT data sources(photometry data + TTL inputs to sync behavior/stimuli)

Assumptions:
        2-color input: 470 nm (calcium dependent signal) and 405 nm (calcium independent control)
        Requires additional standalone functions (provided, credit: Tom Davidson): tdt2mat, controlFit, deltaFF

***********************************************
2) Batch analysis of individual data

What you can expect this script to do:
        Identify timestamps for key behavioral events that are separated in time to avoid re-sampling the same data twice
        Generate PSTH graphs for these events for each subject - batch processing for exploratory visualization of       
            individual variability

Key output data stored in structure arrays: 
        Photometry_PSTH 
        Behavior_ts
        Behavior PSTH (some programs)

 What you need to customize:
        Input/output data filenames
        Specify relationship between TTL inputs and behavioral events
        Specify time window for analysis

Assumptions:
        Raw data have been previously extracted and normalized as in step (1)
   
Requires two additional functions (provided): 
        processPhotDataRow_normDat (to align photometry data across trials to generate PSTH, credit: Tom Davidson)
        opacity (to generate fills for shaded error bars for graphs)
       
***********************************************
3) Generate group average PSTH graphs

What you can expect this script to do:
        Generate group average graphs of behavioral data and photometry responses as PSTH

What you need to customize:
        Input/output data filenames
        Specify time windows for analysis

Assumptions:
      Raw data have been previously analyzed through steps 1 and 2

Requires one additional function (provided): 
       opacity (to generate fills for shaded error bars for graphs)

***********************************************
4) Statistical analysis of group graphs

What you can expect this script to do:
        Compare pre/post event behavior or neural activity using paired t-tests
        Print summary results in excel spreadsheet

Key output data is stored in a structure array "Statistics" with 9 fields:
        Name - Description of comparison
        Comparitor1 - name of data1
        Comparitor2 - name of data2
        Hypothesis - 1 if null hypothesis can be rejected at the 5% level, 0 if not
        pValue - obvious
        Mean1 - group mean for the first comparitor
        Mean2 - group mean for the second comparitor
        Values1 - comparitor 1 individual data used to generate group average
        Values2 - comparitor 2 individual data used to generate group average

 What you need to customize:
       Input/output data filenames
       Specify analysis window

Assumptions:
        Raw data have been previously analyzed through steps 1 and 2

Feel free to contact me with questions - I'm happy to clarify anything in the code. [elizabeth.steinberg at gmail]
