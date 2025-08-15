% SUMMARY - MATLAB
%   For each ROI, model the relationship between PET SUVR and a chosen ASL metric (CBF or ATT),
%   and produce region-wise scatter plots with linear fit and 95% CI, stratified by amyloid status.
%
%   Input: multi-sheet Excel (one sheet per ROI or similar layout)
%   Uses columns: Region, ROI_SUVR, Perfusion_Mean (CBF), Arrival_Mean (ATT),
%                 Amyloid_Status ('Y'/'N'), Age, Voxel_Count, Tracer, QRISK (continuous)

clear all;
close all;
clc;

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

% ---- Renames ----
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

% ---- Set metrics & plotting colours ----
regions = unique(string(data.Region));
colours = struct('CBF', [0, 0.6, 0], 'ATT', [0, 0.4, 1]); % green & blue

metrics      = {'Perfusion_Mean', 'Arrival_Mean'};
metricLabels = {'CBF', 'ATT'};

for m = 1:length(metrics)
    metric = metrics{m};
    label  = metricLabels{m};

    for r = 1:length(regions)
        region  = regions{r};
        subdata = data(strcmp(string(data.Region), region), :);

        % ---- Linear model: ROI_SUVR ~ Metric ----
        mdl = fitlm(subdata, sprintf('ROI_SUVR ~ %s + Tracer', metric), ... % Add covaiates here Age + Tracer + QRISK
                    'CategoricalVars','Tracer');

        coeff = mdl.Coefficients.Estimate(2); % 2nd row corresponds to the 'metric'
        pval  = mdl.Coefficients.pValue(2);
        adjR2 = mdl.Rsquared.Adjusted;
        dof   = mdl.DFE;

        % ---- Scatter plot ----
        figure('Name', sprintf('%s - %s', region, label), 'Color', 'w');
        hold on;

        isPositive = false(height(subdata),1);
        if ismember('Amyloid_Status', subdata.Properties.VariableNames)
            isPositive = strcmp(string(subdata.Amyloid_Status), 'Y');
        end

        scatter(subdata{isPositive,   metric}, subdata{isPositive,   'ROI_SUVR'}, 50, ...
            'V', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours.(label));
        scatter(subdata{~isPositive,  metric}, subdata{~isPositive,  'ROI_SUVR'}, 25, ...
            'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours.(label));

        % ---- Prediction & 95% CI (fix covariates at subgroup means; tracer at mode) ----
        xmin = min(subdata{:, metric}, [], 'omitnan');
        xmax = max(subdata{:, metric}, [], 'omitnan');
        if ~isfinite(xmin) || ~isfinite(xmax) || xmin==xmax
            xmin = min(subdata{:, metric}); xmax = max(subdata{:, metric});
        end
        xvals = linspace(xmin, xmax, 100)';

        % Choose a reference tracer (most frequent in this region)
        if ~isempty(subdata.Tracer)
            refTracer = mode(subdata.Tracer);
        else
            refTracer = categorical("Flutemetamol");
        end

        preddata = table( ...
            xvals, ...
            repmat(mean(subdata.Age,         'omitnan'), 100, 1), ...
            repmat(refTracer,                              100, 1), ...
            repmat(mean(subdata.Qrisk,       'omitnan'), 100, 1), ...
            'VariableNames', {metric, 'Age', 'Tracer', 'Qrisk'} );

        [ypred, yci] = predict(mdl, preddata);

        % CI shading
        fill([xvals; flipud(xvals)], [yci(:,1); flipud(yci(:,2))], ...
            colours.(label), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        % Regression line
        plot(xvals, ypred, '-', 'Color', colours.(label), 'LineWidth', 2);

        % Labels & title
        if strcmp(label, 'CBF')
            xlabel('CBF (mL/100g/min)')
        elseif strcmp(label, 'ATT')
            xlabel('ATT (s)')
        end
        ylabel('Mean SUVR (Aβ PET)');

        title(sprintf('%s \nβ = %.2f, Adj. R^2 = %.2f, dof = %d, pval = %.4f', ...
            region, coeff, adjR2, dof, pval));

        legend('Amyloid Positive', 'Amyloid Negative', '95% CI', 'Regression Line', ...
            'Location', 'best');
        grid on; hold off;
    end
end
