================================================================================
 NM-MRI Processing Pipeline
================================================================================

OVERVIEW
--------
This pipeline processes neuromelanin-sensitive MRI (NM-MRI) data organized in
BIDS format. It produces a single consolidated CSV with raw ROI signal
intensities and a data-driven reference region value for each subject/session.

ROIs extracted: Substantia Nigra (SN, SNc, SNr), VTA, PBP, Locus Coeruleus (LC)


--------------------------------------------------------------------------------
DIRECTORY STRUCTURE
--------------------------------------------------------------------------------

NM_pipeline/
  run_NM_pipeline.m          Main driver script -- run this
  compute_NM_reference.m     Data-driven reference region function (CV-based)
  README.txt                 This file

  NM_toolbox/                NM_toolbox by Wengler et al. (NeuroImage 2020)
    batch_bids.m             Full pipeline (realign, normalize, smooth, CNR)
    batch_bids_raw.m         Raw intensity extraction only
    batch_calc_NM_MNI.m      CNR calculation
    batch_calc_NM_MNI_raw.m  Raw ROI intensity extraction
    batch_normalize_ANTs*.m  ANTs-based MNI normalization variants
    batch_realign.m          Motion correction
    batch_smooth.m           Gaussian smoothing
    default_par.m            Default pipeline parameters
    par_bids.m               BIDS directory parsing
    ... (additional toolbox scripts)
    Masks/                   ROI masks in MNI space
      SN_mask.nii            Substantia Nigra (full)
      SNc_mask_50.nii        SNc (50% probability)
      SNr_mask_50.nii        SNr (50% probability)
      VTA_mask_50.nii        VTA (50% probability)
      PBP_mask_50.nii        PBP (50% probability)
      LC_mask_50.nii         Locus Coeruleus (50% probability)
      CC_mask.nii            Crus Cerebri (reference, used by CNR)
      Overinclusive_LC_Mask*.nii  LC masks for normalization step
      PT_MNI.nii             Pontine tegmentum
    ANTs_Images/             Brain extraction and registration templates
      BrainExtractionBrain_MNI_T1.nii
      Kirby/, NKI/, OASIS/, NKI10andUnder/

  nm_ref_tool/               Data-driven reference region tool
    BrainExtractionMask_MNI_T1_removeROIs.nii  Brain mask (midbrain excluded)
    grow_min_region_intensity_distance.m        Region-growing algorithm
    minN.m                                      N-dimensional minimum utility
    extract_NM_reference.m                      Old F-test version (archived)

  archive_scripts/           Original 3-step scripts (archived, not used)
    step00_move_T1w_files.m
    step01_run_batch_bids.m
    step02_extract_raw.m

  NM_output/                 Pipeline outputs written here (auto-created)
    NM_pipeline_results.csv  Final consolidated output (one row per subject)
    idv/                     Per-subject intermediate CSVs
    reference_region/        Reference region mask + QC maps + log

  _to_delete/                Duplicate files staged for deletion
    ANTs_Images/             (duplicate of NM_toolbox/ANTs_Images/)
    Masks/                   (duplicate of NM_toolbox/Masks/)
    nm_anova/                (empty, from old pipeline)


--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------

  - MATLAB
  - SPM12           (must be on MATLAB path before running)
  - ANTs            (must be available on the system PATH)
  - NM_toolbox      (included in NM_toolbox/)


--------------------------------------------------------------------------------
HOW TO RUN
--------------------------------------------------------------------------------

1. Open run_NM_pipeline.m in MATLAB.

2. Edit the CONFIGURATION block at the top of the file:

     bids_dir     = '/path/to/your/BIDS/root';
     output_dir   = '/path/to/output/directory';
     toolbox_dir  = '/path/to/NM_toolbox';
     ref_tool_dir = '/path/to/nm_ref_tool';

   All other parameters (structfilter, nmfilters, intensity_pct) have
   sensible defaults and can be left as-is unless your dataset requires
   different file naming conventions.

3. Make sure SPM12 is on your MATLAB path:
     addpath('/path/to/spm12')

4. Run the script:
     run_NM_pipeline

   The script adds all required directories to the path automatically.

Expected runtime: dominated by the ANTs normalization step (~10-30 min per
subject depending on hardware). All other steps are fast.


--------------------------------------------------------------------------------
PIPELINE STEPS (what run_NM_pipeline.m does)
--------------------------------------------------------------------------------

STEP 00 -- Copy T1w into nm_data/
  Finds each subject/session's T1w MPRAGE in the BIDS anat/ folder and copies
  it into the corresponding nm_data/ folder so that NM and T1w images are
  co-located for co-registration. Skips files already present.

STEP 01 -- NM_toolbox batch pipeline
  Runs the full NM_toolbox preprocessing chain for all subjects/sessions:
    - Motion correction (realignment to first NM image)
    - Average across NM runs
    - Spatial normalization to MNI space via ANTs
    - Gaussian smoothing (1mm FWHM default)
  Output: smoothed, spatially normalized NM images (swmeanr* files)

STEP 02 -- Raw ROI intensity extraction
  Applies MNI-space ROI masks to the smoothed NM images and extracts mean
  and SD intensity for each ROI (SN, SNc, SNr, VTA, PBP, LC). Raw intensities
  are preferred over the toolbox's CNR metric. Positive-voxel thresholding
  (>0) is applied per ROI.

STEP 03 -- Data-driven reference region (compute_NM_reference.m)
  Constructs a 700-voxel reference region that is:
    (a) Outside the midbrain and ROIs (enforced by BrainExtractionMask_MNI_T1_removeROIs.nii)
    (b) Low intensity (voxels above the 25th percentile of mean intensity excluded)
    (c) Maximally stable across all subjects (minimum Coefficient of Variation)
  Algorithm:
    1. Load all subjects' smoothed NM images
    2. Apply exclusion mask (removes midbrain, ensures searchable brain only)
    3. Apply intensity percentile threshold
    4. Compute voxelwise CV (SD / mean) across subjects
    5. Seed from the minimum-CV voxel
    6. Grow to 700 voxels (balancing CV value and centroid distance)
    7. Validate: flag any subject with zero/NaN coverage in the mask
  Saves: NM_reference_mask.nii, QC intensity/CV maps, and a log file.

STEP 04 -- Consolidated CSV output
  Combines all results into a single file: NM_output/NM_pipeline_results.csv
  Parses subject_id (sub-XXXX) and session (ses-YY) from BIDS paths.


--------------------------------------------------------------------------------
OUTPUT CSV COLUMNS
--------------------------------------------------------------------------------

  subject_id      Subject identifier (e.g. sub-0001)
  session         Session identifier (e.g. ses-01)
  SN_mean_raw     SN mean intensity, all voxels (no threshold)
  SN_mean         SN mean intensity, positive voxels only
  SNc_mean        SNc mean intensity
  SNr_mean        SNr mean intensity
  VTA_mean        VTA mean intensity
  PBP_mean        PBP mean intensity
  LC_mean         LC mean intensity
  SN_stdev        SN stdev, all voxels
  SN_sd           SN stdev, positive voxels
  SNc_sd          SNc stdev
  SNr_sd          SNr stdev
  VTA_sd          VTA stdev
  PBP_sd          PBP stdev
  LC_sd           LC stdev
  NM_reference    Mean intensity in data-driven reference region


--------------------------------------------------------------------------------
BIDS INPUT STRUCTURE EXPECTED
--------------------------------------------------------------------------------

  bids_dir/
    sub-XXXX/
      ses-YY/
        anat/
          sub-XXXX_ses-YY_T1w.nii.gz      (T1w MPRAGE)
        nm_data/
          sub-XXXX_ses-YY_nm*.nii.gz      (NM-MRI scan)

  Note: The T1w file is copied from anat/ into nm_data/ automatically
  by Step 00.


--------------------------------------------------------------------------------
NOTES
--------------------------------------------------------------------------------

- force_rerun flag in run_NM_pipeline.m: set to true to re-run all steps
  even if outputs already exist. Default is false (skip existing outputs).

- intensity_pct in run_NM_pipeline.m: controls how aggressively the reference
  region is constrained to low-intensity voxels. Default is 25 (25th
  percentile). Lower values = more conservative (quieter voxels only).

- The reference region is computed once across all subjects and is fully
  group-agnostic. It does not require a grouping variable and can be reused
  across analyses with different group definitions.

- If any subjects have zero or NaN values in the reference region mask, they
  are flagged in NM_output/reference_region/NM_reference_log.txt and their
  NM_reference value will be NaN in the output CSV.
================================================================================
