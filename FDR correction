clear all;
close all;
clc;

% --- Load table ---
basedir  = 'D:\Yasmin_Liz\Needed\Correct';
datafile = fullfile(basedir, '1Q2Rregression_covar_Qrisk.xlsx');

T = readtable(datafile,'VariableNamingRule','preserve');

T.p_FDR_BH    = nan(height(T),1);
T.sig_FDR_q05 = false(height(T),1);
q = 0.05;

% --- One family per Metric (CBF, ATT, etc.) ---
[G, metricNames] = findgroups(string(T.Metric));

for g = 1:max(G)
    rows = (G==g);
    p = T.pValue(rows);

    v = isfinite(p) & p>0 & p<=1;
    p_adj = nan(size(p));

    if any(v)
% ---- FDR adjusted p-values ----
        [ps, order] = sort(p(v), 'ascend');
        m = numel(ps);

% --- Raw FDR multipliers ---
        p_FDR = ps .* m ./ (1:m)';   % vector, length m

% --- Enforce monotone non-decreasing adjusted p-values ---
        for i = m-1:-1:1
            p_FDR(i) = min(p_FDR(i), p_FDR(i+1));
        end

% --- Map back to original rows; clip at 1 element-wise ---
        tmp = nan(m,1);
        tmp(order) = min(p_FDR, 1);  % element-wise min with 1
        p_adj(v) = tmp;
    end

% --- Write back per-row ---
    T.p_FDR_BH(rows)    = p_adj;
    T.sig_FDR_q05(rows) = p_adj < q;
end

% Stars
stars = strings(height(T),1);
stars(T.p_FDR_BH < 0.05)  = "*";
stars(T.p_FDR_BH < 0.01)  = "**";
stars(T.p_FDR_BH < 0.001) = "***";
T.sig_label = stars;

writetable(T,'Q2_fdr_tracer.xlsx','WriteMode','overwritesheet');

% Check grouping and uniqueness
[G, metricNames] = findgroups(string(T.Metric));
for g = 1:max(G)
    rows = (G==g);
    p_raw = T.pValue(rows);
    p_adj = T.p_FDR_BH(rows);
    fprintf('\nMetric: %s | m=%d | unique raw p=%d | unique FDR p=%d\n', ...
        metricNames(g), sum(rows), numel(unique(round(p_raw,12))), numel(unique(round(p_adj,12))));
end
