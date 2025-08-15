clear all;
close all;
clc;

basedir = 'D:\Yasmin_Liz\DATA'
table_file   = fullfile(basedir, 'FreeSurfer_LD_YT.xlsx');
demo_file = fullfile(basedir, 'Demographics.xlsx');
out_excel = fullfile(basedir, '1Raw_Mass_Table.xlsx')

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

% --- Read the demographics Excel table ---
demo_tbl = readtable(demo_file, 'ReadVariableNames', true, 'VariableNamingRule','preserve');
demo_tbl.Properties.VariableNames = matlab.lang.makeValidName(demo_tbl.Properties.VariableNames);
disp(demo_tbl.Properties.VariableNames);
if isnumeric(demo_tbl.ID)
    demo_tbl.ID = compose('%03d', demo_tbl.ID)
end
demo_IDs = string(demo_tbl.ID)
demo_ages = demo_tbl.Age_yrs_;
demo_amy_status = string(demo_tbl.Amyloid)

% --- Find subject folder ---
subject_dirs = dir(fullfile(basedir, '*_V1*'));
nSubs = numel(subject_dirs);
roi_tables = struct();

% --- Composite definitions ---
ACC_members = {'rostralanteriorcingulate','caudalanteriorcingulate'};
MTL_members = {'left_hippocampus','right_hippocampus','lefthippocampus','righthippocampus','entorhinal'};
% (listed both underscored and un-underscored hippocampus to be robust)

% Initialise composite summary
summary_data = [];

% --- Participant loop ---
for idx = 1:nSubs 
% --- Define paths ---
    subj_folder = subject_dirs(idx).name;
    subj_dir = fullfile(basedir, subj_folder);
    subj_id = extractBefore(subj_folder, '_V1');
    out_dir = fullfile(subj_dir, 'Real_ROI_Extracted_Map_Wholecereb');
    if ~exist(out_dir, 'dir')
    mkdir(out_dir);
    end


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

    T1_file    = dir(fullfile(subj_dir, '*T1w.nii')); %%%% only verifying in same space, not being used for anything
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

% --- Check dimensions ---
     if isequal(size(atlas_img), size(T1_img), size(amy_img), size(perf_img), size(att_img)) == 0
        fprintf('Dimensions mismatch for %s - skipping .\n', subj_id);
        continue;
     end

% ---- Preallocate maps and storage table ----
    medPerf_map = nan(size(atlas_img));
    meanPerf_map = nan(size(atlas_img));
    medArr_map = nan(size(atlas_img));
    meanArr_map = nan(size(atlas_img));
    medAmy_map = nan(size(atlas_img));
    meanAmy_map = nan(size(atlas_img));
    rows = {};

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
    
% --- Weighted mean ---
    nL1 = nnz(cerebL1_mask); nL2 = nnz(cerebL2_mask);
    nR1 = nnz(cerebR1_mask); nR2 = nnz(cerebR2_mask);
    cereb_mean = (nL1*cerebL1_mean + nL2*cerebL2_mean + nR1*cerebR1_mean + nR2*cerebR2_mean) / (nL1+nL2+nR1+nR2);

    cereb_median = median([cerebL1_mean, cerebR1_mean, cerebL2_mean, cerebR2_mean]);

% --- preallocate variables ---
    
    roi_amy_means = zeros(size(roi_indices));
    roi_perf_means = zeros(size(roi_indices));
    roi_perf_median = zeros(size(roi_indices));
    roi_att_means = zeros(size(roi_indices));
    roi_att_median = zeros(size(roi_indices));
    roi_voxcount = zeros(size(roi_indices));

% --- Loop each label (ensures label exists) ---
    for i = 1:length(roi_indices)
        roi_val = roi_indices(i);
        match_roi_labels = find(roi_indices == roi_val, 1);
        roi_name = matlab.lang.makeValidName(roi_labels{match_roi_labels});
        roi_mask = (atlas_img == roi_val);
        if nnz(roi_mask) <= 10, 
            continue; 
        end % skips missing/small ROIs

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
        voxcount_pet  = numel(roi_amy_values);


% --- Compute mean & median --- 
        perfusion_mean_value(i) = mean(roi_perf_values); 
        perfusion_median(i) = median(roi_perf_values);
        arrival_mean_value(i) = mean(roi_att_values);   
        arrival_median(i) = median(roi_att_values);
        amyloid_mean_value(i) = mean(roi_amy_values);   
        amyloid_median(i) = median(roi_amy_values);
        roi_voxcount(i) = nnz(roi_mask);

        fprintf('Computed averages for subject %s\n', subj_id);

% --- Fill median maps ---
        meanPerf_map(roi_mask) = perfusion_mean_value(i);
        medPerf_map(roi_mask) = perfusion_median(i);
        meanArr_map(roi_mask)  = arrival_mean_value(i);
        medArr_map(roi_mask)  = arrival_median(i);
        meanAmy_map(roi_mask)  = amyloid_mean_value(i);
        medAmy_map(roi_mask)  = amyloid_median(i);

        roi_suvr(i) = amyloid_mean_value(i) / cereb_mean;

        fprintf('Subject %s: ROI %s | Cerebellum mean - %.4f |ROI SUVR = %.4f\n', subj_id, roi_name, cereb_mean, roi_suvr(i));

% --- Store metrics for this subject ---
        new_row = table({subj_id}, age, string(amy_stat), voxcount_t1, voxcount_pet, voxcount_cbf, voxcount_att, ...
            amyloid_mean_value(i), amyloid_median(i), perfusion_mean_value(i), perfusion_median(i), arrival_mean_value(i), ...
            arrival_median(i), roi_suvr(i), ...
    'VariableNames', {'ID', 'Age', 'Amyloid Status', 'VoxelCount_T1', 'VoxelCount_PET', 'VoxelCount_CBF', 'VoxelCount_ATT', ...
    'Amyloid Mean', 'Amyloid Median', 'Perfusion Mean', 'Perfusion Median', 'Arrival Mean', 'Arrival Median', 'SUVR'});

        % Append to roi_tables.(roi_name)
        if isfield(roi_tables, roi_name)
            roi_tables.(roi_name) = [roi_tables.(roi_name); new_row];
        else
            roi_tables.(roi_name) = new_row;
        end

    end

% --- Save mean/median maps ---
    out = perf_vol;  out.fname = fullfile(out_dir,'Mean_Perfusion_Map.nii'); spm_write_vol(out, meanPerf_map);
    out = perf_vol;  out.fname = fullfile(out_dir,'Median_Perfusion_Map.nii'); spm_write_vol(out, medPerf_map);
    out = att_vol;   out.fname = fullfile(out_dir,'Mean_Arrival_Map.nii');   spm_write_vol(out, meanArr_map);
    out = att_vol;   out.fname = fullfile(out_dir,'Median_Arrival_Map.nii');   spm_write_vol(out, medArr_map);
    out = amy_vol;   out.fname = fullfile(out_dir,'Mean_Amyloid_Map.nii');  spm_write_vol(out, meanAmy_map);
    out = amy_vol;   out.fname = fullfile(out_dir,'Median_Amyloid_Map.nii');  spm_write_vol(out, medAmy_map);

% --- Save normalised mean/median maps ---
    meanAmy_map_norm = meanAmy_map / cereb_mean;
    medAmy_map_norm = medAmy_map / cereb_median;
    out = amy_vol;   out.fname = fullfile(out_dir,'Mean_Amyloid_Map_Normalised.nii');  spm_write_vol(out, meanAmy_map_norm);
    out = amy_vol;   out.fname = fullfile(out_dir,'Median_Amyloid_Map_Normalised.nii');  spm_write_vol(out, medAmy_map_norm);

end

 roi_list = fieldnames(roi_tables);

for i = 1:numel(roi_list)
    raw_name = roi_list{i};
    sheetname = raw_name;

    % Replace invalid characters with underscore
    sheetname = regexprep(sheetname, '[:\\/*?[\]]', '_');

    % Truncate to 31 characters
    if length(sheetname) > 31
        sheetname = sheetname(1:31);
    end

    this_tbl = roi_tables.(raw_name);
    this_tbl = sortrows(this_tbl, 'SUVR', 'descend');
    writetable(this_tbl, out_excel, 'Sheet', sheetname);
end
