clear all;
close all;
clc;

basedir = 'D:\Yasmin_Liz\DATA'
table_file   = fullfile(basedir, 'FreeSurfer_LD_YT.xlsx');
demo_file = fullfile(basedir, 'Demographics.xlsx');
out_excel = fullfile(basedir, 'Patient_Metrics1.xlsx')

% --- Load the FreeSurfer ROI info ---
roi_table = readtable(table_file, 'VariableNamingRule', 'preserve');
roi_indices   = roi_table.index;
roi_labels    = roi_table.('structure/region name');
roi_inclusion = roi_table.('Inclusion/Exclusion');

% --- Define target cerebellar ROI names ---
cerebL1_index = 8;
cerebL2_index = 7;
cerebR1_index = 47;
cerebR2_index = 46;

include_mask = strcmpi(strtrim(roi_inclusion), 'Inclusion');
included_indices = roi_indices(include_mask);

% --- Read the demographics Excel table ---
 demo_tbl = readtable(demo_file, 'ReadVariableNames', true, 'VariableNamingRule','preserve');
 demo_tbl.Properties.VariableNames = matlab.lang.makeValidName(demo_tbl.Properties.VariableNames);
 disp(demo_tbl.Properties.VariableNames);
if isnumeric(demo_tbl.ID)
    demo_tbl.ID = compose('%03d', demo_tbl.ID)
end
demo_IDs = string(demo_tbl.ID)
demo_ages = demo_tbl.Age_yrs_;
demo_amy_status = string(demo_tbl.Amyloid) %%%%% NEW

% --- Find subject folder ---
subject_dirs = dir(fullfile(basedir, '*_V1*'));
nSubs = numel(subject_dirs);
rows = {};

% --- Participant loop ---
for idx = 1:nSubs
% --- Define paths ---
    subj_folder = subject_dirs(idx).name;
    subj_dir = fullfile(basedir, subj_folder);
    subj_id = extractBefore(subj_folder, '_V1');

% --- Matching participant demograpghics --- 
    match_idx = find(demo_IDs == subj_id);
    age = demo_ages(match_idx);
    amy_stat = demo_amy_status(match_idx);

    fprintf('\nProcessing subject %s...\n', subj_id);

% --- Read files ---
    atlas_file = dir(fullfile(subj_dir, 'aparc.DKTatlas+aseg-in-rawavg_mgz2nii.nii'));
    atlas_vol  = spm_vol(fullfile(subj_dir, atlas_file(1).name));
    atlas_img  = spm_read_vols(atlas_vol);

    amy_file   = dir(fullfile(subj_dir, 'PVC_PET.nii'));
        if isempty(amy_file)
           fprintf('No amyloid image found for %s - skipped.\n', subj_folder);
           continue;
        else 
        amy_vol    = spm_vol(fullfile(subj_dir, 'PVC_PET.nii'));
        amy_img    = spm_read_vols(amy_vol);
        end

    T1_file    = dir(fullfile(subj_dir, '*T1w.nii'));
    if isempty(T1_file)
           fprintf('No T1 image found for %s - skipped.\n', subj_folder);
           T1_img = [];
        else 
    T1_vol     = spm_vol(fullfile(subj_dir, T1_file(1).name));
    T1_img     = spm_read_vols(T1_vol);
    end

    perf_file  = dir(fullfile(subj_dir, 'perfusion_calib_padded_to_T1.nii'));
    if isempty(perf_file)
           fprintf('No perf image found for %s - skipped.\n', subj_folder);
           perf_img = [];
        else 
    perf_vol   = spm_vol(fullfile(subj_dir, perf_file(1).name));
    perf_img   = spm_read_vols(perf_vol);
    end

    att_file   = dir(fullfile(subj_dir, 'arrival_padded_to_T1.nii'));
    if isempty(att_file)
           fprintf('No att image found for %s - skipped.\n', subj_folder);
           att_img = [];
        else 
    att_vol    = spm_vol(fullfile(subj_dir, att_file(1).name));
    att_img    = spm_read_vols(att_vol);
    end

    fprintf('Read images for subject %s\n', subj_id);

% --- Calculate mean amyloid in cerebellum ---
    cerebL1_mask = (atlas_img == cerebL1_index);
    cerebL2_mask = (atlas_img == cerebL2_index);
    cerebR1_mask = (atlas_img == cerebR1_index);
    cerebR2_mask = (atlas_img == cerebR2_index);
    if ~any(cerebL1_mask(:)) || ~any(cerebR1_mask(:))
        error('Cerebellum ROIs not found for subject %s', subj_id);
    end

    cerebL1_values = amy_img(cerebL1_mask);
    cerebL2_values = amy_img(cerebL2_mask);
    cerebR1_values = amy_img(cerebR1_mask);
    cerebR2_values = amy_img(cerebR2_mask);

    if any(isnan(cerebL1_values)) || any(isnan(cerebR1_values))
        error('NaN detected in PET values for cerebellum in subject %s', subj_id);
    end
    
    cerebL1_mean = mean(cerebL1_values);
    cerebL2_mean = mean(cerebL2_values);
    cerebR1_mean = mean(cerebR1_values);
    cerebR2_mean = mean(cerebR2_values);
    
    nL1 = nnz(cerebL1_mask); nL2 = nnz(cerebL2_mask);
    nR1 = nnz(cerebR1_mask); nR2 = nnz(cerebR2_mask);
    cereb_mean = (nL1*cerebL1_mean + nL2*cerebL2_mean + nR1*cerebR1_mean + nR2*cerebR2_mean) / (nL1+nL2+nR1+nR2);

% --- preallocate variables ---
    
    roi_amy_means = zeros(size(included_indices));
    roi_perf_means = zeros(size(included_indices));
    roi_att_means = zeros(size(included_indices));
    roi_voxcount = zeros(size(included_indices)); 

     for r = 1:length(included_indices)
        roi_idx = included_indices(r);
        match_roi_labels = find(roi_indices == roi_idx, 1); %%%%%%%
        roi_name = matlab.lang.makeValidName(roi_labels{match_roi_labels});
    
        % --- Find mask for this ROI in atlas ---
        roi_mask = (atlas_img == roi_idx);
    
% --- Sanity check: mask must not be empty ---
        if ~any(roi_mask(:))
            error('ROI index %d not found in atlas for subject %s', roi_idx, subj_id);
        end
    
% --- Extract values in ROI mask ---
        roi_amy_values = amy_img(roi_mask);
        roi_perf_values = perf_img(roi_mask);
        roi_att_values = att_img(roi_mask);

% --- Omit zeros for ATT/CBF only ---
        roi_perf_values = roi_perf_values(roi_perf_values ~= 0);
        roi_att_values  = roi_att_values(roi_att_values  ~= 0);
    
% --- Check for unexpected NaNs in image values ---
        if any(isnan(roi_amy_values))
            error('NaN detected in PET values for ROI index %d in subject %s', roi_idx, subj_id);
        end
        if any(isnan(roi_perf_values))
            error('NaN detected in Perf values for ROI index %d in subject %s', roi_idx, subj_id);
        end
        if any(isnan(roi_att_values))
            error('NaN detected in ATT values for ROI index %d in subject %s', roi_idx, subj_id);
        end
    
% --- T1 voxel count (atlas mask size) ---
        voxcount_t1   = nnz(roi_mask); % full mask size from T1/atlas
% --- for weighted composites later ---
        voxcount_cbf = numel(roi_perf_values);
        voxcount_att  = numel(roi_att_values);

% --- Compute and store values ---
        roi_amy_means(r) = mean(roi_amy_values); 
        roi_perf_means(r) = mean(roi_perf_values);
        roi_perf_median(r) = median(roi_perf_values)
        roi_att_means(r) = mean(roi_att_values);
        roi_att_median(r) = median(roi_att_values);
        roi_voxcount(r) = nnz(roi_mask); 

        roi_suvr(r) = roi_amy_means(r) / cereb_mean;

        fprintf('Subject %s: ROI %s | Cerebellum mean - %.4f |ROI SUVR = %.4f\n', subj_id, roi_name, cereb_mean, roi_suvr(r));

        rows_cells{r} = {roi_name, age, string(amy_stat), voxcount_t1, voxcount_cbf, voxcount_att, ...
        roi_amy_means(r), roi_perf_means(r), roi_perf_median(r), roi_att_means(r),  roi_att_median(r), roi_suvr(r)};
    end

    tbl = cell2table(vertcat(rows_cells{:}), ...
     'VariableNames', {'ROI','Age','Amyloid_status', 'VoxelCount_T1', 'VoxelCount_CBF','VoxelCount_ATT', ...
     'Amyloid_Mean', 'Perfusion_Mean','Perfusion_Median', 'Arrival_Mean','Arrival_Median', 'ROI_SUVR'});

    % --- Save summary ---
    excel_file = fullfile(basedir, 'Patient_Metrics1.xlsx');
    writetable(tbl, excel_file, 'Sheet', subj_id);
    
    fprintf('Finished subject %s\n', subj_id);

end