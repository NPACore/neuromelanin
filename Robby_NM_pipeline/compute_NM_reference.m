function [nm_reference, ref_mask] = compute_NM_reference(PAR, mask_path, output_dir, intensity_pct)
% COMPUTE_NM_REFERENCE  Data-driven, group-agnostic NM reference region.
%
% Replaces the F-test/ANOVA approach in extract_NM_reference.m.
% Instead of requiring a grouping variable, this function finds the most
% stable, low-intensity region across ALL subjects using the Coefficient
% of Variation (CV = SD / mean). The search space is constrained by:
%   (1) a brain mask that excludes midbrain ROIs (SN, VTA, LC, etc.)
%   (2) an intensity percentile threshold to keep only quiet voxels
%
% INPUTS
%   PAR            - parameter struct from run_NM_pipeline (has nmdirs,
%                    subjects, nsubs, nnm, nmprefs)
%   mask_path      - full path to BrainExtractionMask_MNI_T1_removeROIs.nii
%   output_dir     - directory to save reference mask and diagnostics
%   intensity_pct  - intensity percentile ceiling (default: 25).
%                    Voxels with mean intensity ABOVE this percentile
%                    across subjects are excluded from seed selection.
%
% OUTPUTS
%   nm_reference   - (nsubs x 1) mean intensity within reference region,
%                    one value per subject. NaN if subject image missing
%                    or reference region invalid for that subject.
%   ref_mask       - (logical 3D) the 700-voxel reference region mask
%
% DEPENDENCIES
%   SPM12, grow_min_region_intensity_distance.m, minN.m

if nargin < 4 || isempty(intensity_pct)
    intensity_pct = 25;
end

fprintf('\n--- compute_NM_reference ---\n');
fprintf('Mask          : %s\n', mask_path);
fprintf('Output dir    : %s\n', output_dir);
fprintf('Intensity pct : %.0f\n', intensity_pct);

%% ---- Load exclusion mask ------------------------------------------------
% Mask = 1 inside searchable brain, 0 in midbrain ROIs and outside brain
if ~exist(mask_path, 'file')
    error('Mask not found: %s', mask_path);
end
mask_vol = spm_vol(mask_path);
mask_img = logical(spm_read_vols(mask_vol));
sz = size(mask_img);
fprintf('Exclusion mask loaded. Searchable voxels: %d\n', nnz(mask_img));

%% ---- Load all subjects smoothed NM images --------------------------------
% Uses swmeanr* files -- same images used for ROI extraction in batch_calc_NM_MNI_raw
nm = 1; % reference region uses first (typically only) NM filter
fprintf('Loading smoothed NM images for %d subjects...\n', PAR.nsubs);

all_imgs    = NaN([sz, PAR.nsubs]);   % 4D stack: [x y z subject]
valid_subs  = true(PAR.nsubs, 1);

for sb = 1:PAR.nsubs
    NM_file = spm_select('FPList', PAR.nmdirs{sb, nm}, ...
                          ['^swmeanr' PAR.nmprefs '.*\.(img|nii)$']);
    if isempty(strtrim(NM_file))
        fprintf('  WARNING: No swmeanr* file for %s -- skipping\n', PAR.subjects{sb});
        valid_subs(sb) = false;
        continue;
    end
    NM_vol = spm_vol(NM_file(1,:));
    NM_img = spm_read_vols(NM_vol);
    all_imgs(:,:,:,sb) = NM_img;
    fprintf('  Loaded: %s\n', strtrim(NM_file));
end

n_valid = sum(valid_subs);
fprintf('%d / %d subjects loaded successfully.\n', n_valid, PAR.nsubs);

if n_valid < 2
    error('Need at least 2 valid subjects to compute CV. Aborting.');
end

%% ---- Compute common coverage mask across all subjects -------------------
% Only include voxels where EVERY subject has a valid (non-zero, non-NaN)
% value. Subjects are scanned with slight tilts so edge voxels differ;
% restricting to the shared intersection ensures no subject can ever have
% a zero or NaN within the reference region.
valid_coverage = all_imgs > 0 & ~isnan(all_imgs);              % [x y z nsubs] logical
common_mask    = all(valid_coverage(:,:,:,valid_subs), 4);      % true only where ALL loaded subjects valid
fprintf('Common coverage voxels (valid in every subject): %d\n', nnz(common_mask));

% Combine with the anatomical exclusion mask (removes midbrain / outside brain)
search_mask = mask_img & common_mask;
fprintf('Final searchable voxels (exclusion + common coverage): %d\n', nnz(search_mask));

if nnz(search_mask) < 700
    error(['Fewer than 700 voxels are valid in every subject AND within the ' ...
           'exclusion mask. Check slice coverage overlap across subjects.']);
end

%% ---- Apply combined search mask ----------------------------------------
% Restrict all images to the shared valid space before computing statistics
excl = repmat(~search_mask, [1 1 1 PAR.nsubs]);
all_imgs(excl) = NaN;

%% ---- Compute voxelwise mean and CV across subjects ----------------------
fprintf('Computing voxelwise mean and CV...\n');

mean_img = mean(all_imgs, 4, 'omitnan');
std_img  = std(all_imgs,  0, 4, 'omitnan');
cv_img   = std_img ./ (abs(mean_img) + eps);  % CV = SD / mean

% Reinforce: any voxel outside the search mask gets NaN (belt-and-suspenders)
cv_img(~search_mask) = NaN;

%% ---- Apply intensity threshold ------------------------------------------
% Keep only voxels in the LOWER intensity_pct of mean signal across the
% search mask. This enforces that the reference region is genuinely low-signal.
masked_means = mean_img(search_mask & ~isnan(mean_img));
thresh = prctile(masked_means, intensity_pct);
fprintf('Mean intensity %.0fth-percentile threshold: %.5f\n', intensity_pct, thresh);

cv_img(mean_img > thresh) = NaN;   % exclude high-intensity voxels from seeding
cv_img(~search_mask)      = NaN;   % reinforce combined mask

n_candidates = nnz(~isnan(cv_img));
fprintf('Candidate voxels after intensity filter: %d\n', n_candidates);

if n_candidates < 700
    error(['Fewer than 700 candidate voxels remain after masking and intensity ' ...
           'filtering. Try raising intensity_pct (currently %.0f%%).'], intensity_pct);
end

%% ---- Find minimum-CV seed voxel ----------------------------------------
fprintf('Finding minimum-CV seed voxel...\n');
[min_cv, seed_idx] = minN(cv_img);
fprintf('Seed voxel   : [%d, %d, %d]\n', seed_idx(1), seed_idx(2), seed_idx(3));
fprintf('Seed CV      : %.5f\n', min_cv);
fprintf('Seed mean int: %.5f\n', mean_img(seed_idx(1), seed_idx(2), seed_idx(3)));

%% ---- Grow reference region to 700 voxels --------------------------------
% Replace NaN with Inf so grow function never selects masked-out voxels
fprintf('Growing reference region to 700 voxels...\n');
cv_for_grow = cv_img;
cv_for_grow(isnan(cv_for_grow)) = Inf;

ref_mask = grow_min_region_intensity_distance(cv_for_grow, seed_idx, 700, 0.7);
fprintf('Reference region grown. Final voxel count: %d\n', nnz(ref_mask));

%% ---- Extract reference value per subject --------------------------------
% Because the region was grown entirely within the common coverage mask,
% every valid subject is guaranteed to have full, non-zero coverage.
% Any failure here indicates an unexpected data problem and errors out.
fprintf('Extracting per-subject reference region means...\n');
nm_reference = NaN(PAR.nsubs, 1);

for sb = 1:PAR.nsubs
    if ~valid_subs(sb)
        continue;  % subject had no swmeanr* file; already NaN
    end
    subj_img = all_imgs(:,:,:,sb);
    ref_vals = subj_img(ref_mask);

    if any(ref_vals == 0) || any(isnan(ref_vals)) || isempty(ref_vals)
        % This should not happen if common_mask was built correctly
        error(['Subject %s has zero or NaN values in the reference region ' ...
               'despite common coverage masking. Check image integrity.'], ...
               PAR.subjects{sb});
    end

    nm_reference(sb) = mean(ref_vals);
end

fprintf('Reference extraction complete. %d / %d subjects have valid NM_reference.\n', ...
        sum(~isnan(nm_reference)), PAR.nsubs);

%% ---- Save reference mask and diagnostics --------------------------------
mkdir(output_dir);

% Save binary reference region mask
ref_mask_path = fullfile(output_dir, 'NM_reference_mask.nii');
ref_vol          = mask_vol;
ref_vol.fname    = ref_mask_path;
ref_vol.dt       = [spm_type('float32'), spm_platform('bigend')];
spm_write_vol(ref_vol, double(ref_mask));
fprintf('Reference mask saved : %s\n', ref_mask_path);

% Save mean intensity and CV maps for QC
mean_path = fullfile(output_dir, 'NM_ref_mean_intensity.nii');
cv_path   = fullfile(output_dir, 'NM_ref_CV.nii');

mean_vol       = mask_vol; mean_vol.fname = mean_path;
mean_vol.dt    = [spm_type('float32'), spm_platform('bigend')];
spm_write_vol(mean_vol, mean_img);

cv_vol         = mask_vol; cv_vol.fname = cv_path;
cv_vol.dt      = [spm_type('float32'), spm_platform('bigend')];
spm_write_vol(cv_vol, cv_img);

fprintf('QC maps saved        : %s\n', mean_path);
fprintf('                     : %s\n', cv_path);

% Save a log of key parameters
log_path = fullfile(output_dir, 'NM_reference_log.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'NM Reference Region Log\n');
fprintf(fid, '=======================\n');
fprintf(fid, 'Date/time        : %s\n', datestr(now));
fprintf(fid, 'Subjects         : %d total, %d with swmeanr* files\n', PAR.nsubs, n_valid);
fprintf(fid, 'Common coverage  : %d voxels valid in every subject\n', nnz(common_mask));
fprintf(fid, 'Search space     : %d voxels (coverage + exclusion mask)\n', nnz(search_mask));
fprintf(fid, 'Exclusion mask   : %s\n', mask_path);
fprintf(fid, 'Intensity pct    : %.0f\n', intensity_pct);
fprintf(fid, 'Intensity thresh : %.5f\n', thresh);
fprintf(fid, 'Candidate voxels : %d\n', n_candidates);
fprintf(fid, 'Seed voxel       : [%d, %d, %d]\n', seed_idx(1), seed_idx(2), seed_idx(3));
fprintf(fid, 'Seed CV          : %.5f\n', min_cv);
fprintf(fid, 'Seed mean int    : %.5f\n', mean_img(seed_idx(1), seed_idx(2), seed_idx(3)));
fprintf(fid, 'Final voxel count: %d\n', nnz(ref_mask));
fprintf(fid, 'Valid NM_ref     : %d / %d subjects\n', sum(~isnan(nm_reference)), PAR.nsubs);
fclose(fid);
fprintf('Log saved            : %s\n', log_path);

fprintf('--- compute_NM_reference complete ---\n\n');
end
