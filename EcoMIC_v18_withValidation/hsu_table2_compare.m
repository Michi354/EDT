function out = hsu_table2_compare(T)
% Recreates your Hsu-style GA-band summary (72–96h, non-TH).
% Prints a compact console line per GA band and returns a struct.

need = {'CO_L_min','BSA','HR','GA','PNA','cooling'};
assert(all(ismember(need, T.Properties.VariableNames)), ...
    'hsu_table2_compare: missing columns: %s', ...
    strjoin(setdiff(need, T.Properties.VariableNames), ', '));

CO_Lmin = T.CO_L_min;          % L/min
BSA     = T.BSA;               % m^2
CI      = CO_Lmin ./ BSA;      % L/min/m^2
HR      = T.HR;                % bpm

% Stroke volume
SV_mL   = (CO_Lmin*1000) ./ max(HR, eps);

mask_hsu = (T.cooling==0) & (T.PNA>=3) & (T.PNA<5);

ga_edges   = [0 28 30 32 34 36 38 41];
ga_labels  = {'≤28','29–30','31–32','33–34','35–36','37–38','39–41'};
[g,~]      = discretize(T.GA, ga_edges);

fmt   = @(m,s) sprintf('%.2f ± %.2f', m, s);
fmtHR = @(m,s) sprintf('%.0f ± %.1f', m, s);
fmtSV = @(m,s) sprintf('%.2f ± %.2f', m, s);

n    = zeros(numel(ga_labels),1);
co_s = strings(numel(ga_labels),1);
ci_s = strings(numel(ga_labels),1);
hr_s = strings(numel(ga_labels),1);
sv_s = strings(numel(ga_labels),1);

for k=1:numel(ga_labels)
    idx = mask_hsu & (g==k);
    n(k) = sum(idx);
    if n(k)==0
        co_s(k) = '--'; ci_s(k)='--'; hr_s(k)='--'; sv_s(k)='--';
    else
        co_s(k) = fmt(mean(CO_Lmin(idx),'omitnan'), std(CO_Lmin(idx),'omitnan'));
        ci_s(k) = fmt(mean(CI(idx),'omitnan'),      std(CI(idx),'omitnan'));
        hr_s(k) = fmtHR(mean(HR(idx),'omitnan'),    std(HR(idx),'omitnan'));
        sv_s(k) = fmtSV(mean(SV_mL(idx),'omitnan'), std(SV_mL(idx),'omitnan'));
    end
end

% Hsu 2016 Table 2 reference (mean ± SD)
hsu_CO    = [0.23 0.29 0.35 0.35 0.43 0.47 0.53];
hsu_CO_sd = [0.03 0.06 0.07 0.07 0.08 0.10 0.14];
hsu_CI    = [2.31 2.45 2.69 2.54 2.64 2.49 2.60];
hsu_CI_sd = [0.26 0.24 0.36 0.32 0.33 0.39 0.54];
hsu_HR    = [149 145 142 142 136 131 126];
hsu_HR_sd = [11.4 8.7 8.9 15.7 13.5 11.9 12.1];
hsu_SV    = [1.56 1.99 2.53 2.49 3.22 3.67 4.24];
hsu_SV_sd = [0.28 0.44 0.47 0.60 0.70 0.72 0.89];

hsu_co_s = arrayfun(@(m,s) sprintf('%.2f ± %.2f',m,s), hsu_CO,  hsu_CO_sd, 'uni',0);
hsu_ci_s = arrayfun(@(m,s) sprintf('%.2f ± %.2f',m,s), hsu_CI,  hsu_CI_sd, 'uni',0);
hsu_hr_s = arrayfun(@(m,s) sprintf('%.0f ± %.1f',m,s), hsu_HR,  hsu_HR_sd, 'uni',0);
hsu_sv_s = arrayfun(@(m,s) sprintf('%.2f ± %.2f',m,s), hsu_SV,  hsu_SV_sd, 'uni',0);

fprintf('\nHSU Table 2 comparison (72–96 h, non-TH), by GA group:\n');
for k=1:numel(ga_labels)
    fprintf('  GA %s (n=%d): CO %s vs Hsu %s | CI %s vs Hsu %s | HR %s vs Hsu %s | SV %s vs Hsu %s\n', ...
        ga_labels{k}, n(k), co_s(k), hsu_co_s{k}, ci_s(k), hsu_ci_s{k}, hr_s(k), hsu_hr_s{k}, sv_s(k), hsu_sv_s{k});
end

out = struct('labels',{ga_labels}, 'n',n, 'CO',co_s, 'CI',ci_s, 'HR',hr_s, 'SV',sv_s);
end