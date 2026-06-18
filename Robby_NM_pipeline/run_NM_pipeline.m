%% run_NM_pipeline.m
% Unified NM-MRI pipeline driver
%
% Runs all processing steps in sequence and produces a single consolidated
% CSV with one row per subject/session.
%
% STEPS
%   00  Copy T1w files from anat/ into nm_data/ for co-registration
%   01  Run NM_toolbox batch pipeline (realign, avg, ANTs normalize, smooth)
%   02  Extract raw ROI intensities (SN, SNc, SNr, VTA, PBP, LC)
%   03  Compute data-driven reference region via CV (no grouping variable)
%   04  Consolidate all results into one CSV
%
% OUTPUT CSV COLUMNS
%   subject_id    - parsed from BIDS path (e.g. sub-0001)
%   session       - parsed from BIDS path (e.g. ses-01)
%   SN_mean_raw   - SN mean intensity (all voxels, no threshold)
%   SN_mean       - SN mean intensity (positive voxels only)
%   SNc_mean      - SNc mean intensity
%   SNr_mean      - SNr mean intensity
%   VTA_mean      - VTA mean intensity
%   PBP_mean      - PBP mean intensity
%   LC_mean       - LC mean intensity
%   SN_stdev      - SN stdev (all voxels)
%   SN_sd         - SN stdev (positive voxels)
%   SNc_sd        - SNc stdev
%   SNr_sd        - SNr stdev
%   VTA_sd        - VTA stdev
%   PBP_sd        - PBP stdev
%   LC_sd         - LC stdev
%   SN_CNR_raw    - SN CNR (all voxels, no threshold); CNR = (ROI - CC_mode)/CC_mode * 100
%   SN_CNR        - SN CNR (positive voxels only)
%   SNc_CNR       - SNc CNR
%   SNr_CNR       - SNr CNR
%   VTA_CNR       - VTA CNR
%   PBP_CNR       - PBP CNR
%   LC_CNR        - LC CNR
%   SN_CNR_stdev  - SN CNR stdev (all voxels)
%   CC_mode       - Histogram mode of CC intensity (value used to normalise CNR)
%   CC_sd         - CC intensity standard deviation
%   CC_mean       - Arithmetic mean of CC intensity
%   NM_reference  - Mean intensity in data-driven CV-based reference region
%
% DEPENDENCIES
%   MATLAB, SPM12 (on path), ANTs (on system PATH)
%   NM_toolbox functions      (toolbox_dir)
%   compute_NM_reference.m    (same directory as this script)
%   grow_min_region_intensity_distance.m, minN.m  (ref_tool_dir)

%% =========================================================================
%% CONFIGURATION -- edit these paths for your dataset
%% =========================================================================

bids_dir     = '/data/D3/Nifti';               % BIDS root (contains sub-*/ses-*/anat/ etc.)
output_dir   = '/data/NEUROMELANIN/NM_output'; % where all outputs are written
toolbox_dir  = '/data/NEUROMELANIN/NM_toolbox'; % NM_toolbox directory (batch_bids, default_par, etc.)
ref_tool_dir = '/data/NEUROMELANIN/nm_ref_tool'; % nm_ref_tool directory (grow_min, minN, mask)

% NM_toolbox file matching
structfilter = 'T1w.nii.gz';          % T1w filename pattern
nmfilters    = {'nm*.nii.gz'};        % NM filename pattern(s)

% Reference region settings
mask_path     = fullfile(ref_tool_dir, 'BrainExtractionMask_MNI_T1_removeROIs.nii');
intensity_pct = 25;   % only voxels below this percentile of mean intensity
                      % are eligible as reference region (keeps quiet voxels)

% Overwrite behaviour
% Set to true to re-run steps even if outputs already exist
force_rerun = false;

%% =========================================================================
%% DEPENDENCY CHECK
%% =========================================================================

fprintf('\n========================================\n');
fprintf(' NM-MRI Pipeline\n');
fprintf(' %s\n', datestr(now));
fprintf('========================================\n\n');

if isempty(which('spm'))
    error('SPM12 not found on MATLAB path. Add SPM12 before running.');
end

if ~exist(bids_dir, 'dir')
    error('BIDS directory not found: %s', bids_dir);
end

if ~exist(toolbox_dir, 'dir')
    error('NM_toolbox directory not found: %s', toolbox_dir);
end

if ~exist(ref_tool_dir, 'dir')
    error('nm_ref_tool directory not found: %s', ref_tool_dir);
end

if ~exist(mask_path, 'file')
    error('Exclusion mask not found: %s', mask_path);
end

% Add all required directories to MATLAB path:
%   pipeline_dir  -- so compute_NM_reference.m is findable
%   toolbox_dir   -- NM_toolbox functions (batch_bids, default_par, etc.)
%   ref_tool_dir  -- grow_min_region_intensity_distance.m, minN.m
pipeline_dir = fileparts(mfilename('fullpath')); % directory this script lives in
addpath(pipeline_dir);
addpath(toolbox_dir);
addpath(ref_tool_dir);

% Switch into toolbox dir so default_par can locate itself via mfilename
cd(toolbox_dir);

%% =========================================================================
%% STEP 00 -- Copy T1w files into nm_data/ folders
%% =========================================================================

fprintf('--- STEP 00: Copying T1w files into nm_data/ ---\n');

anat_list = dir(fullfile(bids_dir, '*', '*', 'anat', '*T1w.nii*'));

if isempty(anat_list)
    % Also try without session level (sub-XX/anat/)
    anat_list = dir(fullfile(bids_dir, '*', 'anat', '*T1w.nii*'));
end

n_copied  = 0;
n_skipped = 0;

for idx = 1:length(anat_list)
    anatpath     = fullfile(anat_list(idx).folder, anat_list(idx).name);
    anat_newpath = strrep(anatpath, '/anat/', '/nm_data/');

    if ~exist(anat_newpath, 'file') || force_rerun
        % Ensure nm_data directory exists
        nm_data_dir = fileparts(anat_newpath);
        if ~exist(nm_data_dir, 'dir')
            mkdir(nm_data_dir);
        end
        copyfile(anatpath, anat_newpath);
        fprintf('  Copied : %s\n', anat_list(idx).name);
        n_copied = n_copied + 1;
    else
        n_skipped = n_skipped + 1;
    end
end

fprintf('Step 00 complete: %d copied, %d already present.\n\n', n_copied, n_skipped);

%% =========================================================================
%% BUILD PAR STRUCT (shared by Steps 01-03)
%% =========================================================================

PAR_update.root         = absdir(output_dir);
PAR_update.structfilter = structfilter;
PAR_update.nmfilters    = nmfilters;

PAR = default_par(PAR_update);
PAR = par_bids(PAR, absdir(bids_dir));
PAR = batch_bids2root(PAR);

fprintf('Subjects found: %d\n', PAR.nsubs);
for sb = 1:PAR.nsubs
    fprintf('  %s\n', PAR.subjects{sb});
end
fprintf('\n');

%% =========================================================================
%% STEP 01 -- Run NM_toolbox batch pipeline
%% =========================================================================

fprintf('--- STEP 01: NM_toolbox batch pipeline ---\n');

% Realign NM images (motion correction within subject)
batch_realign(PAR);

% Average motion-corrected NM images
batch_calc_avg(PAR);

% Spatial normalization to MNI via ANTs (slow step)
fprintf('Starting ANTs normalization (this takes a while)...\n');
mni_start = tic();
batch_normalize_ANTs_optRedo(PAR);
fprintf('ANTs normalization complete. Time: %.1f min\n', toc(mni_start)/60);

% Gaussian smoothing of spatially normalized images
batch_smooth(PAR);

fprintf('Step 01 complete.\n\n');

%% =========================================================================
%% STEP 02a -- Extract CNR values (CC-normalised, uses CC_mask.nii)
%% =========================================================================

fprintf('--- STEP 02a: CNR extraction (CC reference) ---\n');

cnr_results = batch_calc_NM_MNI(PAR);
% cnr_results cell array: row 1 = headers, rows 2:end = one per subject
% Columns:
%  1  Subject              2  NM Dir
%  3  SN CNR Raw (mean)    4  SN CNR (mean)    5  SNc CNR (mean)
%  6  SNr CNR (mean)       7  VTA CNR (mean)   8  PBP CNR (mean)   9  LC CNR (mean)
%  10 SN CNR (stdev)
%  11 CC Intensity (mode)   <- hist_mode, used for CNR normalisation
%  12 CC Intensity (stdev)
%  13 CC Intensity (mean)   <- arithmetic mean

fprintf('Step 02a complete.\n\n');

%% =========================================================================
%% STEP 02b -- Extract raw ROI intensities
%% =========================================================================

fprintf('--- STEP 02b: Raw ROI intensity extraction ---\n');

roi_results = batch_calc_NM_MNI_raw(PAR);
% roi_results cell array: row 1 = headers, rows 2:end = one per subject
% Columns:
%  1  Subject           2  NM Dir
%  3  SN_mean_raw       4  SN_mean      5  SNc_mean
%  6  SNr_mean          7  VTA_mean     8  PBP_mean     9  LC_mean
%  10 SN_stdev(all)     11 SN_sd        12 SNc_sd
%  13 SNr_sd            14 VTA_sd       15 PBP_sd       16 LC_sd

fprintf('Step 02b complete.\n\n');

%% =========================================================================
%% STEP 03 -- Compute data-driven reference region
%% =========================================================================

fprintf('--- STEP 03: Data-driven reference region (CV-based) ---\n');

ref_output_dir = fullfile(output_dir, 'reference_region');
[nm_reference, ref_mask] = compute_NM_reference(PAR, mask_path, ref_output_dir, intensity_pct);

fprintf('Step 03 complete.\n\n');

%% =========================================================================
%% STEP 04 -- Build consolidated CSV
%% =========================================================================

fprintf('--- STEP 04: Building consolidated CSV ---\n');

% --- Parse subject_id and session from PAR.subjects ---
% PAR.subjects entries are in BIDS format: 'sub-XXXX/ses-YY' or 'sub-XXXX'
subject_id = cell(PAR.nsubs, 1);
session    = cell(PAR.nsubs, 1);

for sb = 1:PAR.nsubs
    raw = PAR.subjects{sb};
    % Extract sub-* token
    sub_tok = regexp(raw, '(sub-[^/]+)', 'tokens', 'once');
    if ~isempty(sub_tok)
        subject_id{sb} = sub_tok{1};
    else
        subject_id{sb} = raw;  % fallback: use as-is
    end
    % Extract ses-* token if present
    ses_tok = regexp(raw, '(ses-[^/]+)', 'tokens', 'once');
    if ~isempty(ses_tok)
        session{sb} = ses_tok{1};
    else
        session{sb} = '';
    end
end

% --- Pull raw intensity values from roi_results ---
n_rows = PAR.nsubs;

SN_mean_raw = cell2mat(roi_results(2:n_rows+1, 3));
SN_mean     = cell2mat(roi_results(2:n_rows+1, 4));
SNc_mean    = cell2mat(roi_results(2:n_rows+1, 5));
SNr_mean    = cell2mat(roi_results(2:n_rows+1, 6));
VTA_mean    = cell2mat(roi_results(2:n_rows+1, 7));
PBP_mean    = cell2mat(roi_results(2:n_rows+1, 8));
LC_mean     = cell2mat(roi_results(2:n_rows+1, 9));
SN_stdev    = cell2mat(roi_results(2:n_rows+1, 10));
SN_sd       = cell2mat(roi_results(2:n_rows+1, 11));
SNc_sd      = cell2mat(roi_results(2:n_rows+1, 12));
SNr_sd      = cell2mat(roi_results(2:n_rows+1, 13));
VTA_sd      = cell2mat(roi_results(2:n_rows+1, 14));
PBP_sd      = cell2mat(roi_results(2:n_rows+1, 15));
LC_sd       = cell2mat(roi_results(2:n_rows+1, 16));

% --- Pull CNR values from cnr_results ---
SN_CNR_raw   = cell2mat(cnr_results(2:n_rows+1, 3));
SN_CNR       = cell2mat(cnr_results(2:n_rows+1, 4));
SNc_CNR      = cell2mat(cnr_results(2:n_rows+1, 5));
SNr_CNR      = cell2mat(cnr_results(2:n_rows+1, 6));
VTA_CNR      = cell2mat(cnr_results(2:n_rows+1, 7));
PBP_CNR      = cell2mat(cnr_results(2:n_rows+1, 8));
LC_CNR       = cell2mat(cnr_results(2:n_rows+1, 9));
SN_CNR_stdev = cell2mat(cnr_results(2:n_rows+1, 10));
CC_mode      = cell2mat(cnr_results(2:n_rows+1, 11)); % hist_mode of CC intensity (used for CNR)
CC_sd        = cell2mat(cnr_results(2:n_rows+1, 12));
CC_mean      = cell2mat(cnr_results(2:n_rows+1, 13)); % arithmetic mean of CC intensity

% --- Assemble table ---
T = table( ...
    subject_id, session, ...
    SN_mean_raw, SN_mean, SNc_mean, SNr_mean, VTA_mean, PBP_mean, LC_mean, ...
    SN_stdev, SN_sd, SNc_sd, SNr_sd, VTA_sd, PBP_sd, LC_sd, ...
    SN_CNR_raw, SN_CNR, SNc_CNR, SNr_CNR, VTA_CNR, PBP_CNR, LC_CNR, ...
    SN_CNR_stdev, CC_mode, CC_sd, CC_mean, ...
    nm_reference, ...
    'VariableNames', { ...
        'subject_id', 'session', ...
        'SN_mean_raw', 'SN_mean', 'SNc_mean', 'SNr_mean', 'VTA_mean', 'PBP_mean', 'LC_mean', ...
        'SN_stdev', 'SN_sd', 'SNc_sd', 'SNr_sd', 'VTA_sd', 'PBP_sd', 'LC_sd', ...
        'SN_CNR_raw', 'SN_CNR', 'SNc_CNR', 'SNr_CNR', 'VTA_CNR', 'PBP_CNR', 'LC_CNR', ...
        'SN_CNR_stdev', 'CC_mode', 'CC_sd', 'CC_mean', ...
        'NM_reference' ...
    });

% --- Write CSV ---
csv_path = fullfile(output_dir, 'NM_pipeline_results.csv');
writetable(T, csv_path);

fprintf('Consolidated CSV saved: %s\n', csv_path);
fprintf('Rows: %d  |  Columns: %d\n', height(T), width(T));
fprintf('\n========================================\n');
fprintf(' Pipeline complete.\n');
fprintf(' %s\n', datestr(now));
fprintf('========================================\n\n');
