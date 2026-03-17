function out = walther_bw_strat(T, texOutFile)
% Week-1, non-TH, 500 g BW bands; prints console summary and writes LaTeX.

need = {'CO_L_min','BW','GA','PNA','cooling'};
assert(all(ismember(need, T.Properties.VariableNames)), ...
    'walther_bw_strat: missing columns: %s', ...
    strjoin(setdiff(need, T.Properties.VariableNames), ', '));

CO_Lmin = T.CO_L_min;
BW_kg   = T.BW/1000;
CO_mlkg = 1000*CO_Lmin./BW_kg;
CO_ml   = CO_Lmin*1000;
PNA_d   = T.PNA;

mask = (T.cooling==0) & (PNA_d <= 7);
T1      = T(mask,:);
CO_mlkg1= CO_mlkg(mask);
CO_ml1  = CO_ml(mask);
BW1     = T1.BW; GA1 = T1.GA; PNA1 = T1.PNA;

bw_edges = 750:500:4750;
bw_edges = bw_edges(bw_edges <= (max(BW1)+500));
bw_labels = strings(numel(bw_edges)-1,1);
for k=1:numel(bw_labels), bw_labels(k) = sprintf('%d–%d', bw_edges(k), bw_edges(k+1)-1); end
g = discretize(BW1, bw_edges);

meanSD = @(x) sprintf('%.0f ± %.0f', mean(x,'omitnan'), std(x,'omitnan'));
medIQR = @(x) sprintf('%.0f [%.0f–%.0f]', median(x,'omitnan'), prctile(x,25), prctile(x,75));

nBin = numel(bw_labels);
n   = zeros(nBin,1);
bw_mean  = strings(nBin,1); ga_mean  = strings(nBin,1); age_mean = strings(nBin,1);
co_ml_s  = strings(nBin,1); co_mlkg_s= strings(nBin,1);
co_ml_med= strings(nBin,1); co_mlkg_med = strings(nBin,1);
pct_walther = nan(nBin,1);

for k=1:nBin
    idx = (g==k);
    n(k) = sum(idx);
    if n(k)>0
        bw_mean(k)      = meanSD(BW1(idx));
        ga_mean(k)      = meanSD(GA1(idx));
        age_mean(k)     = meanSD(PNA1(idx));
        co_ml_s(k)      = meanSD(CO_ml1(idx));
        co_mlkg_s(k)    = meanSD(CO_mlkg1(idx));
        co_ml_med(k)    = medIQR(CO_ml1(idx));
        co_mlkg_med(k)  = medIQR(CO_mlkg1(idx));
        pct_walther(k)  = 100*mean(CO_mlkg1(idx)>=170 & CO_mlkg1(idx)<=327);
    else
        bw_mean(k)='--'; ga_mean(k)='--'; age_mean(k)='--'; co_ml_s(k)='--';
        co_mlkg_s(k)='--'; co_ml_med(k)='--'; co_mlkg_med(k)='--';
    end
end

fprintf('\nWALTHER-like (week-1, non-TH) by 500 g BW bins:\n');
for k=1:nBin
    if n(k)>0
        fprintf('  %s g: n=%d | CO/kg %s (mean±SD), %s (median[IQR]) | Walther 200–324: %.1f%%\n', ...
            bw_labels(k), n(k), co_mlkg_s(k), co_mlkg_med(k), pct_walther(k));
    end
end

% Write LaTeX (optional)
if nargin>=2 && ~isempty(texOutFile)
    fid = fopen(texOutFile,'w');
    assert(fid>0, 'Cannot open %s for writing.', texOutFile);
    fprintf(fid,'\\begin{table}[ht]\n\\footnotesize\n\\centering\n');
    fprintf(fid,'\\begin{tabular}{lrrrrrrr}\n\\hline\n');
    fprintf(fid,'\\textbf{BW group (g)} & \\textbf{n} & \\textbf{BW (g)} & \\textbf{GA (wk)} & \\textbf{Age (d)} & \\textbf{CO (mL/min)} & \\textbf{CO (mL/min/kg)} & \\textbf{%% in 200--324}\\\\\n');
    fprintf(fid,'\\hline\n');
    for k=1:nBin
        if n(k)>0
            fprintf(fid,'%s & %d & %s & %s & %s & %s & %s & %.1f\\\\\n', ...
                bw_labels(k), n(k), bw_mean(k), ga_mean(k), age_mean(k), co_ml_s(k), co_mlkg_s(k), pct_walther(k));
        end
    end
    fprintf(fid,'\\hline\n\\end{tabular}\n');
    fprintf(fid,'\\vspace{0.1cm}\n');
    fprintf(fid,'\\caption{Walther-like stratification: week-1, non-TH cohort summarized by 500~g birth-weight bands.}\n');
    fprintf(fid,'\\label{tab:walther_like}\n\\end{table}\n');
    fclose(fid);
    fprintf('Wrote LaTeX table to: %s\n', fullfile(pwd,texOutFile));
end

out = table(bw_labels, n, bw_mean, ga_mean, age_mean, co_ml_s, co_mlkg_s, co_ml_med, co_mlkg_med, pct_walther, ...
            'VariableNames', {'BWband','n','BW','GA','Age','CO_mLmin','CO_mLmin_kg','CO_median','COkg_median','Pct_200_324'});
end

