%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writes Excel file calculating SUVRs (voxel-weighted)
%
% Reads patient demographics for age, amyloid status and ID matching
% Uses FreeSurfer atlas for ROI definition
% Calculates regional PET means, a voxel-count–weighted cerebellar mean,
% and a voxel-count–weighted target composite mean, then outputs Composite SUVR
% Used for manual Welches t-test performed in MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

basedir     = 'D:\Yasmin_Liz\DATA';
table_file  = fullfile(basedir, 'FreeSurfer_LD_YT.xlsx');
demo_file   = fullfile(basedir, 'Demographics.xlsx');
out_excel   = fullfile(basedir, '1GroupSUVR_Valid.xlsx');

% --- Load the FreeSurfer ROI info ---
roi_table      = readtable(table_file, 'VariableNamingRule', 'preserve');
roi_indices    = roi_table.index;
roi_labels     = roi_table.('structure/region name');
roi_inclusion  = roi_table.('Inclusion/Exclusion');

% --- Define target cerebellar ROI indices (reference region) ---
cerebL1_index = 8;
cerebL2_index = 7;
cerebR1_index = 47;
cerebR2_index = 46;

% --- Which ROIs to include in the TARGET composite (from your sheet) ---
include_mask     = strcmpi(strtrim(roi_inclusion), 'Inclusion');
included_indices = roi_indices(include_mask);

% --- Read demographics ---
demo_tbl = readtable(demo_file, 'ReadVariableNames', true, 'VariableNamingRule','preserve');
demo_tbl.Properties.VariableNames = matlab.lang.makeValidName(demo_tbl.Properties.VariableNames);
disp(demo_tbl.Properties.VariableNames);
if isnumeric(demo_tbl.ID)
    demo_tbl.ID = compose('%03d', demo_tbl.ID);
end
demo_IDs       = string(demo_tbl.ID);
demo_ages      = demo_tbl.Age_yrs_;
demo_amy_stat  = string(demo_tbl.Amyloid);

% --- Subjects ---
subject_dirs  = dir(fullfile(basedir, '*_V1*'));
nSubs         = numel(subject_dirs);
summary_data  = {};  % {ID, Age, Amyloid Status, Composite_SUVR}

for idx = 1:nSubs
    % --- Paths / IDs ---
    subj_folder = subject_dirs(idx).name;
    subj_dir    = fullfile(basedir, subj_folder);
    subj_id     = extractBefore(subj_folder, '_V1');
    fprintf('\nProcessing subject %s...\n', subj_id);

    % --- Atlas (label image in T1 space) ---
    atlas_file = dir(fullfile(subj_dir, 'aparc.DKTatlas+aseg-in-rawavg_mgz2nii.nii'));
    if isempty(atlas_file)
        fprintf('No atlas for %s - skipped.\n', subj_folder);
        continue;
    end
    atlas_vol = spm_vol(fullfile(subj_dir, atlas_file(1).name));
    atlas_img = spm_read_vols(atlas_vol);

    % --- PET image (no PVC version here) ---
    amy_file = dir(fullfile(subj_dir, 'PVC_PET.nii'));
    if isempty(amy_file)
        fprintf('No amyloid image found for %s - skipped.\n', subj_folder);
        continue;
    else
        amy_vol = spm_vol(fullfile(subj_dir, ['PVC_PET.nii']));
        amy_img = spm_read_vols(amy_vol);
    end

    % Included regions
    % We'll compute per-ROI mean and voxel count (PET voxels), then weight.

    roi_means    = zeros(size(included_indices));
    roi_voxcount = zeros(size(included_indices));

    for r = 1:length(included_indices)
        roi_idx = included_indices(r);
        roi_mask = (atlas_img == roi_idx);

        if ~any(roi_mask(:))
            error('ROI index %d not found in atlas for subject %s', roi_idx, subj_id);
        end

        % PET values within ROI (keep zeros; do NOT filter)
        roi_values = amy_img(roi_mask);

        % NaN check only (do not omit NaNs silently)
        if any(isnan(roi_values))
            error('NaN detected in PET values for ROI index %d in subject %s', roi_idx, subj_id);
        end

        roi_means(r)    = mean(roi_values);      % mean over all PET voxels in ROI
        roi_voxcount(r) = numel(roi_values);     % PET voxel count (for weighting)
    end

    % Cerebellum reference
    cerebL1_mask = (atlas_img == cerebL1_index);
    cerebL2_mask = (atlas_img == cerebL2_index);
    cerebR1_mask = (atlas_img == cerebR1_index);
    cerebR2_mask = (atlas_img == cerebR2_index);

    if ~any(cerebL1_mask(:)) || ~any(cerebR1_mask(:))
        error('Cerebellum ROI not found for subject %s', subj_id);
    end

    cerebL1_values = amy_img(cerebL1_mask);
    cerebL2_values = amy_img(cerebL2_mask);
    cerebR1_values = amy_img(cerebR1_mask);
    cerebR2_values = amy_img(cerebR2_mask);

    if any(isnan(cerebL1_values)) || any(isnan(cerebR1_values)) || ...
       any(isnan(cerebL2_values)) || any(isnan(cerebR2_values))
        error('NaN detected in PET values for cerebellum in subject %s', subj_id);
    end

    cerebL1_mean = mean(cerebL1_values);
    cerebL2_mean = mean(cerebL2_values);
    cerebR1_mean = mean(cerebR1_values);
    cerebR2_mean = mean(cerebR2_values);

    nL1 = numel(cerebL1_values);
    nL2 = numel(cerebL2_values);
    nR1 = numel(cerebR1_values);
    nR2 = numel(cerebR2_values);

    % Weighted mean across cerebellar segments (voxel-count weighting)
    cereb_mean = (nL1*cerebL1_mean + nL2*cerebL2_mean + nR1*cerebR1_mean + nR2*cerebR2_mean) ...
                 / (nL1 + nL2 + nR1 + nR2);

% composite (weighted)
    % Weighted target mean over all included ROIs:
    total_vox            = sum(roi_voxcount);
    target_mean_weighted = sum(roi_means .* roi_voxcount) / total_vox;

    % Normalise by cerebellar mean to get composite SUVR:
    composite_suvr = target_mean_weighted / cereb_mean;

    fprintf('Subject %s: Cerebellum mean = %.4f | Composite SUVR (weighted) = %.4f\n', ...
            subj_id, cereb_mean, composite_suvr);

    % --- Demographics lookup ---
    match_idx   = find(demo_IDs == subj_id);
    if isempty(match_idx)
        age        = NaN;
        amy_status = string(missing);
    else
        age        = demo_ages(match_idx);
        amy_status = demo_amy_stat(match_idx);
    end

    % --- Collect summary row ---
    summary_data(end+1,:) = {subj_id, age, amy_status, composite_suvr};
end

% --- Final table & write ---
summary_tbl = cell2table(summary_data, 'VariableNames', {'ID','Age','Amyloid Status','Composite_SUVR'});
summary_tbl = sortrows(summary_tbl, 'Composite_SUVR', 'descend');
writetable(summary_tbl, out_excel);

fprintf('Saved composite SUVR table to %s\n', out_excel);

