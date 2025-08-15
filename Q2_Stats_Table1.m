clear all;
close all;
clc;

% Load data
basedir  = 'D:\Yasmin_Liz\Needed\Correct';
datafile = fullfile(basedir, '1Raw_Roi_Table.xlsx');
sheetNames = sheetnames(datafile);

allData = table();
for k = 1:numel(sheetNames)
    T = readtable(datafile, 'Sheet', sheetNames{k}, 'VariableNamingRule','preserve');

    % Add Region column if missing
    if ~ismember('Region', T.Properties.VariableNames)
        T.Region = repmat(string(sheetNames{k}), height(T), 1);
    end

    allData = [allData; T]; 
end

% Now assign to data for consistency
data = allData;

% ---- Light-touch renames to match code expectations ----
if ismember('Amyloid Status', data.Properties.VariableNames)
    data.Properties.VariableNames{strcmp(data.Properties.VariableNames,'Amyloid Status')} = 'Amyloid_Status';
end

if ismember('Qrisk', data.Properties.VariableNames)
    data.Properties.VariableNames{strcmp(data.Properties.VariableNames,'Qrisk')} = 'Qrisk';
end

% ---- Tracer ----
if ~iscategorical(data.Tracer)
data.Tracer = categorical(data.Tracer);
end

regions = unique(data.Region);
metrics = {'Perfusion_Mean', 'Arrival_Mean'};
metricLabels = {'CBF', 'ATT'};

% Preallocate output
results = struct('Region', {}, 'Metric', {}, 'Beta', {}, 'CI_Lower', {}, 'CI_Upper', {}, 'pValue', {}, 'Adjusted_R2', {},...
    'DoF', {}, 'Covariate', {}, 'Cov_Beta', {}, 'Cov_CI_Lower', {}, 'Cov_CI_Upper', {}, 'Cov_pValue', {});

covName = 'Tracer';   % or Age or Qrisk

% Loop over metrics and regions
for m = 1:length(metrics)
    metric = metrics{m};
    label = metricLabels{m};
   
    for r = 1:length(regions)
        region = regions{r};
        subdata = data(strcmp(data.Region, region), :);
       
        % Fit regression
        mdl = fitlm(subdata, sprintf('ROI_SUVR ~ %s + %s', metric, covName));

 % Extract coefficient + stats
        coeff = mdl.Coefficients.Estimate(2);
        CI = coefCI(mdl); % 95% CI
        ci_lower = CI(2,1);
        ci_upper = CI(2,2);
        pval = mdl.Coefficients.pValue(2);
        adjR2 = mdl.Rsquared.Adjusted;
        dof = mdl.DFE;

        if strcmp(covName, 'Tracer')
            % finds the non-reference tracer term
            cov_idx = find(startsWith(mdl.CoefficientNames, 'Tracer_'));
        else
            % Age/ Qrisk match by name
            cov_idx = find(strcmp(mdl.CoefficientNames, covName));
        end

 % --- covariate stats ---
        cov_beta  = mdl.Coefficients.Estimate(cov_idx);
        cov_p     = mdl.Coefficients.pValue(cov_idx);
        cov_ci    = CI(cov_idx,:); 

        % Add to struct
        results(end+1) = struct( ...
            'Region', region, ...
            'Metric', label, ...
            'Beta', coeff, ...
            'CI_Lower', ci_lower, ...
            'CI_Upper', ci_upper, ...
            'pValue', pval, ...
            'Adjusted_R2', adjR2, ...
            'DoF', dof, ...
            'Covariate', covName, ...
            'Cov_Beta', cov_beta, ...
            'Cov_CI_Lower', cov_ci(1), ...
            'Cov_CI_Upper', cov_ci(2), ...
            'Cov_pValue', cov_p ...
        );
    end
end

% Convert to table
resultTable = struct2table(results);

% Write to Excel
writetable(resultTable, '1Q2Rregression_covar_Tracer.xlsx');
disp('Saved: 1Q2Rregression_covar_Tracer.xlsx');