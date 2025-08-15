% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writes Excel file summarising ROI metrics for all participants
%
% Reads patients demographics file for age, amyloid status and ID matching
% ROIs are define by 'Inclusion/Exclusion' label from read in excel sheet
% Uses Freesurfer atlas for ROI definition, T1 for voxel count, PET, and ASL metrics (mean and median)
% Calculates Mean Normalised PET values
% Output is a single sheet per (L+R combined) region with all participants metrics for that region
% CHANGES:
%   - Robust normalisation of region names: strip ctx/dkt prefix, lh/rh, and left/right (any separator)
%   - Union masks:
%         * ACC = rostral anterior cingulate ∪ caudal anterior cingulate (both hemispheres)
%         * Medial_Temporal = hippocampus (L+R) ∪ entorhinal
%   - Do NOT output separate hippocampus or entorhinal sheets; only the Medial_Temporal composite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

basedir = 'D:\Yasmin_Liz\DATA';
table_file   = fullfile(basedir, 'FreeSurfer_LD_YT.xlsx');
demo_file = fullfile(basedir, 'Demographics.xlsx');
out_excel = fullfile(basedir, '1Raw_Roi_Table.xlsx');

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
    demo_tbl.ID = compose('%03d', demo_tbl.ID);
end
demo_IDs = string(demo_tbl.ID);
demo_ages = demo_tbl.Age_yrs_;
demo_amy_status = string(demo_tbl.Amyloid);

% --- Find subject folder ---
subject_dirs = dir(fullfile(basedir, '*_V1*'));
nSubs = numel(subject_dirs);
roi_tables = struct();

% --- Participant loop ---
for idx = 1:nSubs
    % --- Define paths ---
    subj_folder = subject_dirs(idx).name;
    subj_dir = fullfile(basedir, subj_folder);
    subj_id = extractBefore(subj_folder, '_V1');

    % --- Matching participant demographics ---
    match_idx = find(demo_IDs == subj_id);
    if isempty(match_idx)
        fprintf('No demographics match for %s, skipping.\n', subj_id);
        continue
    end
    age = demo_ages(match_idx);
    amy_stat = demo_amy_status(match_idx);

    fprintf('\nProcessing subject %s...\n', subj_id);

    % --- Read in files ---
    atlas_file = dir(fullfile(subj_dir, 'aparc.DKTatlas+aseg-in-rawavg_mgz2nii.nii'));
    if isempty(atlas_file)
        fprintf('No atlas for %s - skipped.\n', subj_folder);
        continue
    end
    atlas_vol  = spm_vol(fullfile(subj_dir, atlas_file(1).name));
    atlas_img  = spm_read_vols(atlas_vol);

    amy_file   = dir(fullfile(subj_dir, 'PVC_PET.nii'));
    if isempty(amy_file)
        fprintf('No amyloid image found for %s - skipped.\n', subj_folder);
        continue
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
        error('Cerebellum ROI not found for subject %s', subj_id);
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

    % --- Accumulator for pooled values (per region) ---
    acc = struct();  % acc.(region_key).amy / .cbf / .att / .vox_t1 / .disp

    % --- Init union masks for composites ---
    hip_mask = false(size(atlas_img));      % hippocampus (L+R)
    ent_mask = false(size(atlas_img));      % entorhinal
    acc_mask = false(size(atlas_img));      % rostral+caudal ACC

    % --- included ROI indices and pool values ---
    for r = 1:length(included_indices)
        roi_idx = included_indices(r);
        match_roi_labels = find(roi_indices == roi_idx, 1);
        roi_label_raw = roi_labels{match_roi_labels};
        roi_name_full = matlab.lang.makeValidName(roi_label_raw);

        % --- Mask for this ROI in atlas ---
        roi_mask = (atlas_img == roi_idx);
        if ~any(roi_mask(:))
            error('ROI index %d not found in atlas for subject %s', roi_idx, subj_id);
        end

        % ======================= Normalise names  =======================
        % Start lowercased raw, remove optional leading "ctx dkt", then hemispheric tags:
        %   - textual: left/right (+ optional _ or -)
        %   - short:   lh/rh (+ optional _ or -)
        base_raw_lower = lower(strtrim(roi_label_raw));
        base_no_ctx = regexprep(base_raw_lower, '^(ctx[\s\-_]*dkt[\s\-_]*)', '');
        base_no_hemi = regexprep(base_no_ctx, '^(lh|rh|left|right)[\s\-_]*', '');
        base_disp = regexprep(strtrim(base_no_hemi), '^(?i)(lh|rh|left|right)[\s\-_]*', ''); % display-friendly
        base_key  = matlab.lang.makeValidName(strtrim(base_no_hemi));

        % --- Build union masks for composites (use normalised names) ---
        if strcmp(base_key, 'hippocampus')
            hip_mask = hip_mask | roi_mask;
        end
        if strcmp(base_key, 'entorhinal')
            ent_mask = ent_mask | roi_mask;
        end
        is_acc_piece = startsWith(base_key,'caudalanteriorcingulate') || startsWith(base_key,'rostralanteriorcingulate');
        if is_acc_piece
            acc_mask = acc_mask | roi_mask;
        end

        % --- Decide if we add an individual sheet for this region ---
        % Skip rostral/caudal ACC parts (we will write only the composite ACC).
        if is_acc_piece
            continue
        end
        % Skip hippocampus and entorhinal individual sheets (only Medial_Temporal composite wanted).
        if strcmp(base_key,'hippocampus') || strcmp(base_key,'entorhinal')
            continue
        end

        % --- Use the normalised name as the region key ---
        region_key = base_key;
        region_disp = regexprep(base_disp,'[\s\-]+',' ');

        % --- Extract values in ROI mask ---
        roi_amy_values  = amy_img(roi_mask);   % PET: zeros allowed
        roi_perf_values = [];
        roi_att_values  = [];
        if ~isempty(perf_img), roi_perf_values = perf_img(roi_mask); end
        if ~isempty(att_img),  roi_att_values  = att_img(roi_mask);  end

        % --- Omit zeros for ATT/CBF only ---
        roi_perf_values = roi_perf_values(roi_perf_values ~= 0);
        roi_att_values  = roi_att_values( roi_att_values  ~= 0);

        % --- Accumulate pooled data per region ---
        if ~isfield(acc, region_key)
            acc.(region_key).amy    = [];
            acc.(region_key).cbf    = [];
            acc.(region_key).att    = [];
            acc.(region_key).vox_t1 = 0;
            acc.(region_key).disp   = region_disp;
        end
        acc.(region_key).amy    = [acc.(region_key).amy;  roi_amy_values(:)];
        acc.(region_key).cbf    = [acc.(region_key).cbf;  roi_perf_values(:)];
        acc.(region_key).att    = [acc.(region_key).att;  roi_att_values(:)];
        acc.(region_key).vox_t1 = acc.(region_key).vox_t1 + nnz(roi_mask);
    end

    % ------------------------ ADD COMPOSITE SHEETS VIA UNION MASKS ------------------------
    % ACC composite (rostral+caudal anterior cingulate)
    if any(acc_mask(:))
        vals_pet = amy_img(acc_mask);
        if ~isempty(perf_img)
            vals_cbf = perf_img(acc_mask); vals_cbf = vals_cbf(vals_cbf~=0);
        else
            vals_cbf = [];
        end
        if ~isempty(att_img)
            vals_att = att_img(acc_mask);  vals_att = vals_att(vals_att~=0);
        else
            vals_att = [];
        end

        voxcount_t1  = nnz(acc_mask);
        voxcount_pet = numel(vals_pet);
        voxcount_cbf = numel(vals_cbf);
        voxcount_att = numel(vals_att);

        amy_mean = mean(vals_pet);
        perf_mean = mean(vals_cbf);
        perf_med  = median(vals_cbf);
        att_mean  = mean(vals_att);
        att_med   = median(vals_att);

        roi_suvr_val = amy_mean / cereb_mean;

        new_row_acc = table({subj_id}, age, string(amy_stat), ...
            voxcount_t1, voxcount_pet, voxcount_cbf, voxcount_att, ...
            amy_mean, perf_mean, perf_med, att_mean, att_med, roi_suvr_val, ...
            'VariableNames', {'ID','Age','Amyloid Status', ...
            'VoxelCount_T1','VoxelCount_PET','VoxelCount_CBF','VoxelCount_ATT', ...
            'Amyloid Mean','Perfusion_Mean','Perfusion_Median', ...
            'Arrival_Mean','Arrival_Median','ROI_SUVR'});

        if isfield(roi_tables, 'ACC')
            roi_tables.ACC = [roi_tables.ACC; new_row_acc];
        else
            roi_tables.ACC = new_row_acc;
        end
    end

    % Medial_Temporal composite (hippocampus union entorhinal)
    mt_mask = hip_mask | ent_mask;
    if any(mt_mask(:))
        vals_pet = amy_img(mt_mask);
        if ~isempty(perf_img)
            vals_cbf = perf_img(mt_mask); vals_cbf = vals_cbf(vals_cbf~=0);
        else
            vals_cbf = [];
        end
        if ~isempty(att_img)
            vals_att = att_img(mt_mask);  vals_att = vals_att(vals_att~=0);
        else
            vals_att = [];
        end

        voxcount_t1  = nnz(mt_mask);
        voxcount_pet = numel(vals_pet);
        voxcount_cbf = numel(vals_cbf);
        voxcount_att = numel(vals_att);

        amy_mean = mean(vals_pet);
        perf_mean = mean(vals_cbf);
        perf_med  = median(vals_cbf);
        att_mean  = mean(vals_att);
        att_med   = median(vals_att);

        roi_suvr_val = amy_mean / cereb_mean;

        new_row_mt = table({subj_id}, age, string(amy_stat), ...
            voxcount_t1, voxcount_pet, voxcount_cbf, voxcount_att, ...
            amy_mean, perf_mean, perf_med, att_mean, att_med, roi_suvr_val, ...
            'VariableNames', {'ID','Age','Amyloid Status', ...
            'VoxelCount_T1','VoxelCount_PET','VoxelCount_CBF','VoxelCount_ATT', ...
            'Amyloid Mean','Perfusion_Mean','Perfusion_Median', ...
            'Arrival_Mean','Arrival_Median','ROI_SUVR'});

        if isfield(roi_tables, 'Medial_Temporal')
            roi_tables.Medial_Temporal = [roi_tables.Medial_Temporal; new_row_mt];
        else
            roi_tables.Medial_Temporal = new_row_mt;
        end
    end

    % --- After pooling: compute metrics once per kept region and append to roi_tables ---
    base_list = fieldnames(acc);
    for b = 1:numel(base_list)
        bkey = base_list{b};
        amy_vec = acc.(bkey).amy;      % PET (zeros kept)
        cbf_vec = acc.(bkey).cbf;      % CBF (zeros removed)
        att_vec = acc.(bkey).att;      % ATT (zeros removed)

        % counts
        voxcount_t1  = acc.(bkey).vox_t1;
        voxcount_pet = numel(amy_vec);
        voxcount_cbf = numel(cbf_vec);
        voxcount_att = numel(att_vec);

        % metrics (means/medians on pooled vectors)
        amy_mean = mean(amy_vec);
        perf_mean = mean(cbf_vec);
        perf_med  = median(cbf_vec);
        att_mean  = mean(att_vec);
        att_med   = median(att_vec);

        roi_suvr_val = amy_mean / cereb_mean;

        % --- Stores metrics for this subject ---
        new_row = table({subj_id}, age, string(amy_stat), ...
            voxcount_t1, voxcount_pet, voxcount_cbf, voxcount_att, ...
            amy_mean, perf_mean, perf_med, att_mean, att_med, roi_suvr_val, ...
            'VariableNames', {'ID','Age','Amyloid Status', ...
            'VoxelCount_T1','VoxelCount_PET','VoxelCount_CBF','VoxelCount_ATT', ...
            'Amyloid Mean','Perfusion_Mean','Perfusion_Median', ...
            'Arrival_Mean','Arrival_Median','ROI_SUVR'});

        if isfield(roi_tables, bkey)
            roi_tables.(bkey) = [roi_tables.(bkey); new_row];
        else
            roi_tables.(bkey) = new_row;
        end
    end

end

% --- Write out the stored roi_tables in MATLAB to Excel sheet ---
roi_list = fieldnames(roi_tables);
for i = 1:numel(roi_list)
    raw_name = roi_list{i};
    sheetname = raw_name;

    % --- Replace invalid characters with underscores ---
    sheetname = regexprep(sheetname, '[:\\/*?[\]]', '_');

    % Truncate to 31 characters
    if length(sheetname) > 31
        sheetname = sheetname(1:31);
    end

    this_tbl = roi_tables.(raw_name);
    this_tbl = sortrows(this_tbl, 'ROI_SUVR', 'descend');
    writetable(this_tbl, out_excel, 'Sheet', sheetname);
end
