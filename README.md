Multimodal PET–MRI (ASL) — MSc Project

This repository contains MATLAB scripts used in my MSc dissertation on relationships between amyloid PET burden (SUVR) and ASL MRI metrics — cerebral blood flow (CBF) and arterial transit time (ATT) — across FreeSurfer-defined cortical regions. PET values are cerebellum-normalised and analyses include group comparisons and regression models with covariates and interactions.

Data & folder expectations

One folder per subject: <ID>_V1/

Required files per subject (when used by a script):
*_T1w.nii (structural T1)
aparc.DKTatlas+aseg-in-rawavg_mgz2nii.nii (FreeSurfer/DKT atlas in T1 space)
PET_Ave_flipped.nii (pre-averaged PET) and/or rPET_Ave_flipped.nii (PET resliced to T1- used for with no PVC comparison)
PVC_PET.nii (PVC PET; used in ROI extraction script)
perfusion_calib_padded_to_T1.nii (CBF)
arrival_padded_to_T1.nii (ATT)

Support tables in the project root:
FreeSurfer_LD_YT.xlsx (ROI indices, region labels, inclusion flags)
Demographics.xlsx (ID, Age, Amyloid, QRISK_mean_BP)
Cerebellar reference indices (from FreeSurfer_LD_YT.xlsx): 7, 8, 46, 47 (cortex and white matter).

Pipeline at a glance:
Coregistration: PET → T1 (SPM, estimate+write).
ROI extraction & SUVR:
Compute ROI means/medians and cerebellum-normalised SUVRs (per ROI, per subject).
Export per-ROI Excel sheets and summary NIfTI maps.
Compute composite SUVR (mean of included ROI SUVRs) for each subject.

Region-wise aggregation:

Create per-subject, per-region tables (optionally voxel-count–weighted means).

Statistics & visualisation:
Welch’s t-tests (Y vs N).
Linear regressions: unadjusted, covariates, and interactions with Amyloid status or QRISK.
Region-wise scatter plots (SUVR vs ASL metric) with regression line and 95% CI.

Note: OLD files included as analysis pipeline developed differently demanding different outputs and processing/data handling requirements
