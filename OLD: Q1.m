% - MATLAB
% script loads regional PET SUVR data from a specified Excel file and 
% performs region-wise statistical analysis. For each ROI, it fits a linear 
% regression model to assess the relationship between amyloid status and SUVR, 
% extracting beta coefficients, confidence intervals, p-values, and model fit 
% statistics. Outputs are stored in a structured array for further analysis.

clear all;
close all;
clc;

basedir = 'D:\Yasmin_Liz\DATA\Cerebellum_No_Cun_or_medial_temp\Stats_output';
% Load data
datafile = fullfile(basedir, 'Stats_wholecereb.xlsx');
data = readtable('Stats_wholecereb.xlsx');

regions = unique(data.Region);
Group = {'Amyloid_Status'};
GroupNames = data.Group;
isPositive = strcmp(GroupNames, 'Y');
Group1 = data.Group(isY);
Group1 = data.Group(~isY);

% Preallocate output
results = struct('Region', {}, 'pValue', {}, 't_stat', {}, 'CI_Lower', {}, 'CI_Upper', {},  'MeanY', {},'MeanN', {},'Percent_Diff', {});

 for r = 1:length(regions)
        region = regions{r};
        subdata = data(strcmp(data.Region, region), :);
       
        % Fit regression
        [H,P,CI,STATS] = ttest2()]
        mdl = fitlm(subdata, sprintf('Mean_SUVR ~ %s', metric));
       
        % Extract coefficient + stats
        coeff = mdl.Coefficients.Estimate(2);
        CI = coefCI(mdl); % 95% CI
        ci_lower = CI(2,1);
        ci_upper = CI(2,2);
        pval = mdl.Coefficients.pValue(2);
        adjR2 = mdl.Rsquared.Adjusted;
        dof = mdl.DFE;

        % Add to struct
        results(end+1) = struct( ...
            'Region', region, ...
            'Metric', label, ...
            'Beta', coeff, ...
            'CI_Lower', ci_lower, ...
            'CI_Upper', ci_upper, ...
            'pValue', pval, ...
            'Adjusted_R2', adjR2, ...
            'DoF', dof ...
        );

    end

