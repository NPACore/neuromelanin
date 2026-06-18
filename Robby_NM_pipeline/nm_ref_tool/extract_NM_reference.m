function extract_NM_reference(table_fname,group_variable,NM_variable,NM_reference_region_variable,save_dir,mask)
%% -------- DESCRIPTION --------
% Function loads in the table which has 1 row per subject: (1) it has a group
% variable that is a numeric value that indicates which group that subject
% belongs to; (2) paths for each subject's processed NM file. The function
% then conducts an ANOVA and extracts an F-test in spm that tests the null
% hypothesis that any group differs from any other group. It then uses this
% to find the voxel with the lowest difference. It then grows that to 700
% voxels for a reference region. It does this by finding its nearest
% neighbor to the centroid of the current mask and lowest intensity
% neighbor. It then extracts the reference region for each subject and
% saves out the table.

%% -------- INPUTS --------
% table_fname = path to table which has the grouping variable and paths to NM files [string, full path]
% group_variable = name of grouping variable in table [string]
% NM_variable = name of variable that indicates paths of each subjects NM processed file in table [string]
% NM_reference_region_variable = name of reference region in table (this will extract each subject's reference region) [string]
% save_dir = path that indicates where to save anova output [string, full path]
% mask = path to whole brain mask that excludes areas of interest [string, full path]

%% -------- OUTPUTS --------
% Extracts a reference region into the table that differs the least by
% groups on the NM variable. 

%% -------- EXAMPLE --------
% table_fname = 'Y:\Workspaces\Staff\Helmet_Karim\Public\nm\NM_groups.xlsx';
% group_variable = 'Group';
% NM_variable = 'NM_fpath';
% NM_reference_region_variable = 'NM_reference';
% save_dir = 'Y:\Workspaces\Staff\Helmet_Karim\Public\nm\nm_anova\';
% mask = 'Y:\Workspaces\Staff\Helmet_Karim\Public\nm\BrainExtractionMask_MNI_T1_removeROIs.nii';
% extract_NM_reference(table_fname,group_variable,NM_variable,NM_reference_region_variable,save_dir,mask);

%% -------- FUNCTION --------
% load data
data_full = readtable(table_fname);

% remove empty paths
lst = ~cellfun(@isempty,data_full.(NM_variable));
data = data_full(lst,:);

% create directory
mkdir(save_dir);

% run anova
matlabbatch{1}.spm.stats.factorial_design.dir{1} = save_dir;

% create group contrast
group_list = unique(data.(group_variable));
n_groups = length(group_list);
n_contrasts = nchoosek(n_groups, 2);

contrast_matrix = zeros(n_contrasts, n_groups);

% index for each row
row = 1;

% generate all pairs (i, j) where i < j
for i = 1:n_groups-1
    for j = i+1:n_groups
        contrast = zeros(1, n_groups);
        contrast(i) = 1;
        contrast(j) = -1;
        contrast_matrix(row, :) = contrast;
        row = row + 1;
    end
end

% for each unique group, create cells
for idx = 1:length(group_list)
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(idx).scans = data.(NM_variable)(data.(group_variable) == group_list(idx));
end
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = mask;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'group';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = contrast_matrix;
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

save(strcat(save_dir,'anova.mat'),'matlabbatch');
spm_jobman('run',strcat(save_dir,'anova.mat'));

% smooth map
clear matlabbatch;
matlabbatch{1}.spm.spatial.smooth.data{1} = strcat(save_dir,'spmF_0001.nii');
matlabbatch{1}.spm.spatial.smooth.fwhm = [4 4 4];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
save(strcat(save_dir,'smooth.mat'),'matlabbatch');
spm_jobman('run',strcat(save_dir,'smooth.mat'));

% find the minimum value and grow it to 700 voxels
hdr = load_untouch_nii(strcat(save_dir,'sspmF_0001.nii'));
new_hdr = load_untouch_nii(mask);
hdr.img(new_hdr.img == 0) = NaN;
hdr.img(hdr.img < 0.2) = NaN;

[~,idx] = minN(hdr.img);
mask_img = grow_min_region_intensity_distance(hdr.img, idx, 700, 0.7);

new_hdr.img = double(mask_img);
save_untouch_nii(new_hdr,strcat(save_dir,'reference_mask.nii'));

% extract the reference region
data_full.(NM_reference_region_variable) = NaN(size(data_full,1),1);
for subj = 1:size(data_full,1)
    disp(subj);
    if exist(data_full.(NM_variable){subj},'file') & ~isempty(data_full.(NM_variable){subj})
        hdr = load_untouch_nii(data_full.(NM_variable){subj});
        data_full.(NM_reference_region_variable)(subj) = nanmean(hdr.img(mask_img)); %#ok
    end
end

% save
writetable(data_full,table_fname);
end