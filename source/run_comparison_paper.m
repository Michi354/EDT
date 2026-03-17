function run_comparison_paper(models)
%% ========= Paths & cohorts you want to run =================================
files.fn_amicool  = '.\dataverse_files\data as published amicool.xlsx';
files.fn_creacool = '.\dataverse_files\dataset as published creacool.xlsx';
files.fn_nonmem   = '.\dataverse_files\NONMEM.xlsx';

baseOut = './PaperOnOurData';
if ~exist(baseOut,'dir'), mkdir(baseOut); end

% models{flag+1} is used by predict_ours(tab, flag, models)
cohorts = [0 1 2];  % 0=No disease, 1=AKI, 2=Asphyxia+TH

for flag = cohorts
    [tab, label] = load_dataset_by_flag(flag, files);
    outdir = fullfile(baseOut, sprintf('%d_%s', flag, label));
    if ~exist(outdir,'dir'), mkdir(outdir); end

    fprintf('\n================= Cohort %d | %s =================\n', flag, label);

    %% 1) Observed sCr & PAPER mean function
    [obs_colname, y_obs] = getObservedSCr(tab);
    fprintf('Using observed sCr column: %s\n', obs_colname);

    paper_mu = @(GA,PNA) ...
        ( 1 + (0.225 + 0.0102.*GA) .* exp(-0.251 .* (PNA + 0.181.*GA - 8.6).^2) ) ...
      .* ( (53.1 - 1.01.*GA).^3.57 + PNA.^3.57 ) ...
      .* ( 0.966 - 0.0132.*GA ) ...
      ./ ( (53.1 - 1.01.*GA).^3.57 + (1.46 + 0.015.*GA).*PNA.^3.57 );

    %% 2) Predictions: Paper, PaperCal (α,τ), Ours
    y_paper = paper_mu(tab.GA, tab.PNA);

    % robust α,τ on log-scale; internally clamps PNA+τ to >=0 to avoid (neg)^3.57
    [tau, alpha] = calibrate_alpha_tau(tab.GA, tab.PNA, y_obs, paper_mu, -3, 3, 0.05);
    PNA_tau = max(tab.PNA + tau, 0);                     % guard for fractional power
    y_paperCal = exp(alpha) .* paper_mu(tab.GA, PNA_tau);

    y_ours = predict_ours(tab, flag, models);

    %% 3) Common mask (fair head-to-head)
    idx = isfinite(y_obs) & y_obs>0 & isfinite(tab.GA) & isfinite(tab.PNA) ...
        & isfinite(y_paper) & isfinite(y_paperCal) & isfinite(y_ours);
    fprintf('Common-mask n = %d / %d rows kept for metrics.\n', sum(idx), height(tab));
    fprintf('PaperCal calibration: alpha=%.3f (scale=%.3f), tau=%.2f days\n', alpha, exp(alpha), tau);

    %% 4) Metric bundles
    Mp    = metrics_bundle_mask(y_paper,     y_obs, idx);
    Mpc   = metrics_bundle_mask(y_paperCal,  y_obs, idx);
    Mo    = metrics_bundle_mask(y_ours,      y_obs, idx);

    % Full detail -> CSV (not printed)
    MetricsTbl_full = metrics_table_full({'Paper'; 'PaperCal'; 'Ours'}, Mp, Mpc, Mo);
    writetable(MetricsTbl_full, fullfile(outdir,'metrics_extended.csv'));

    % Print only essentials
    PrintTbl = metrics_table_print({'Paper'; 'PaperCal'; 'Ours'}, Mp, Mpc, Mo);
    disp(PrintTbl);

    %% 5) Log–log calibration (intercept≈0, slope≈1 ideal) → print compact
    Cp  = calibration_log_mask(y_paper,     y_obs, idx);
    Cpc = calibration_log_mask(y_paperCal,  y_obs, idx);
    Co  = calibration_log_mask(y_ours,      y_obs, idx);

    CalibTbl_print = table( ...
        {'Paper';'PaperCal';'Ours'}, ...
        [Cp.intercept;  Cpc.intercept;  Co.intercept], ...
        [Cp.slope;      Cpc.slope;      Co.slope], ...
        [Cp.R2;         Cpc.R2;         Co.R2], ...
        'VariableNames',{'Model','Intercept_log','Slope_log','R2'});
    disp(CalibTbl_print);
    writetable(CalibTbl_print, fullfile(outdir,'calibration_log.csv'));

    %% 6) Winner with MAE & RMSE 
    rank_and_print_winner({'Paper','PaperCal','Ours'}, Mp, Mpc, Mo);

    % %% 7) Fairness: scale-only calibration for OURS (print α and write mini-CSV)
    % alpha_ours = mean( log(y_obs(idx)) - log(max(y_ours(idx),eps)) );
    % y_oursCal  = exp(alpha_ours) .* y_ours;
    % M_oCal     = metrics_bundle_mask(y_oursCal, y_obs, idx);
    % fprintf('OURS scale-only calib: alpha_ours=%.3f (scale=%.3f)\n', alpha_ours, exp(alpha_ours));
    % 
    % T_primary = table( ...
    %     {'Paper';'PaperCal';'Ours';'OursCal'}, ...
    %     [Mp.MAPE; Mpc.MAPE; Mo.MAPE; M_oCal.MAPE], ...
    %     'VariableNames', {'Model','MAPE_pct'});
    % writetable(T_primary, fullfile(outdir,'primary_metrics_compare.csv'));
    % fprintf('Saved primary metrics compare to %s\n', fullfile(outdir,'primary_metrics_compare.csv'));

    %% 8) MAPE by observed-y bands (only prints band MAPE for PaperCal vs Ours)
    print_yband_mape(y_obs, y_paperCal, y_ours, idx);

    %% 9) Global scale (Paper→Obs) + empirical PNA-band ratios (saved, not printed)
    [scale_hat, ci_sc] = global_scale_via_logratio(y_obs, tab.GA, tab.PNA, paper_mu);
    fprintf('Global scale (Paper→Obs) = %.3f; 95%% CI [%.3f, %.3f]\n', scale_hat, ci_sc(1), ci_sc(2));

    T_all = empirical_scale_pnaBands(tab, y_obs, paper_mu, [0 1 2 3 4 7 14 21 28], [34 36]);
    writetable(T_all.allGA,    fullfile(outdir,'empirical_scale_pnaBands_allGA.csv'));
    writetable(T_all.gaWindow, fullfile(outdir,'empirical_scale_pnaBands_GA34to36.csv'));
    fprintf('Saved empirical scale tables to %s\n', outdir);

    % Quick ratio plot (saved)
    try
        f = figure('Color','w','Position',[80 80 900 420]); 
        x = 1:height(T_all.allGA);
        errorbar(x, T_all.allGA.Ratio, ...
                 T_all.allGA.Ratio - T_all.allGA.Ratio_lo95, ...
                 T_all.allGA.Ratio_hi95 - T_all.allGA.Ratio, ...
                 'o-','LineWidth',1.5,'MarkerSize',5);
        xticks(x); xticklabels(erase(T_all.allGA.Band, 'PNA ')); xtickangle(0);
        yline(1,'k:'); grid on; box on;
        ylabel('Observed / Paper (median ratio)'); title(sprintf('Empirical scale vs PNA (%s)', label));
        exportgraphics(f, fullfile(outdir,'empirical_scale_vs_PNA_allGA.png'), 'Resolution', 180);
        close(f);
        fprintf('Saved plot: %s\n', fullfile(outdir,'empirical_scale_vs_PNA_allGA.png'));
    catch
        % plotting optional
    end
end

fprintf('\nAll cohorts done.\n');
end

%% ========================== HELPERS (self-contained) ======================

function [tab, label] = load_dataset_by_flag(flag, files)
    switch flag
        case 0, tab = readtable(files.fn_nonmem,  'VariableNamingRule','preserve'); label='NONMEM';
        case 1, tab = readtable(files.fn_creacool,'VariableNamingRule','preserve'); label='Creacool';
        case 2, tab = readtable(files.fn_amicool,'VariableNamingRule','preserve'); label='Amicool';
        otherwise, error('Unknown flag %d',flag);
    end
    % Minimal sanity for GA/PNA
    need = {'GA','PNA'};
    miss = setdiff(need, tab.Properties.VariableNames);
    if ~isempty(miss)
        error('dataset missing required columns: %s', strjoin(miss, ', '));
    end
end

function [tau, alpha] = calibrate_alpha_tau(GA, PNA, y, mu_fun, tau_min, tau_max, tau_step)
% Robust search in log-space: for each tau, closed-form alpha.
    y = double(y);
    m = isfinite(y) & y>0 & isfinite(GA) & isfinite(PNA);
    GA = GA(m); PNA = PNA(m); y = y(m);
    logy = log(y);

    taus = tau_min:tau_step:tau_max;
    bestSSE = Inf; tau=0; alpha=0;
    for t = taus
        mu = mu_fun(GA, PNA + t); mu = max(mu, eps);
        r  = logy - log(mu);
        a  = mean(r);
        SSE = sum( (r - a).^2 );
        if SSE < bestSSE, bestSSE=SSE; tau=t; alpha=a; end
    end
end

function M = metrics_bundle_mask(yhat, y, idx)
% Return a struct with all key metrics (on idx)
    yh = double(yhat(idx)); yo = double(y(idx));
    yh = max(yh, eps);      % guard
    yo = max(yo, eps);

    loge = log(yh) - log(yo);
    abs_e = abs(yh - yo);

    M.N      = numel(yo);
    M.GMAPE  = exp(mean(abs(loge))) - 1;       % geometric MAPE as ratio
    M.GMAPE_pct = 100*M.GMAPE;

    M.MAE    = mean(abs_e);
    M.RMSE   = sqrt(mean((yh - yo).^2));
    M.MedAE  = median(abs_e);

    M.sMAPE  = mean( 2*abs(yh-yo) ./ (abs(yh)+abs(yo)) );
    M.sMAPE_pct = 100*M.sMAPE;

    M.MPE    = mean(yh - yo);
    M.MAPE   = mean(abs((yh-yo)./yo))*100;

    % Calibration-in-the-large (log-normal GMR = geometric mean of yh/yo)
    r   = log(yh ./ yo);
    mu  = mean(r); se = std(r)/sqrt(numel(r));
    M.GMR       = exp(mu);
    ci          = mu + [-1 1]*1.96*se;
    M.GMR_lo95  = exp(ci(1));
    M.GMR_hi95  = exp(ci(2));
end

function T = metrics_table_full(names, Mp, Mpc, Mo)
    row = @(M) [M.N, M.GMAPE_pct, M.MAE, M.RMSE, M.MedAE, ...
                M.sMAPE_pct, M.MPE, M.MAPE, M.GMR, M.GMR_lo95, M.GMR_hi95];
    MetricsMat = [ row(Mp); row(Mpc); row(Mo) ];
    T = array2table(MetricsMat, ...
        'VariableNames', {'N','GMAPE_pct','MAE_mgdl','RMSE_mgdl','MedAE_mgdl', ...
                          'sMAPE_pct','MPE_mgdl','MAPE_pct','GMR','GMR_lo95','GMR_hi95'});
    T.Model = names;
    T = movevars(T, 'Model', 'Before', 1);
end

function T = metrics_table_print(names, Mp, Mpc, Mo)
% Minimal, publication-friendly console print
    row = @(M) [M.N, M.MAE, M.RMSE, M.MAPE, M.GMR, M.GMR_lo95, M.GMR_hi95];
    MetricsMat = [ row(Mp); row(Mpc); row(Mo) ];
    T = array2table(MetricsMat, ...
        'VariableNames', {'N','MAE_mgdl','RMSE_mgdl','MAPE_pct','GMR','GMR_lo95','GMR_hi95'});
    T.Model = names(:);
    T = movevars(T, 'Model', 'Before', 1);
end

function C = calibration_log_mask(yhat, y, idx)
% Regress log(y) on log(yhat): intercept≈0, slope≈1 ideal; report R2
    yh = max(double(yhat(idx)), eps);
    yo = max(double(y(idx)), eps);
    X  = log(yh); Y = log(yo);
    mdl = fitlm(X, Y);
    C.intercept = mdl.Coefficients.Estimate(1);
    C.slope     = mdl.Coefficients.Estimate(2);
    C.R2        = mdl.Rsquared.Ordinary;
end

function rank_and_print_winner(names, Mp, Mpc, Mo)
    MAEs   = [Mp.MAE,    Mpc.MAE,   Mo.MAE];
    RMSEs  = [Mp.RMSE,   Mpc.RMSE,  Mo.RMSE];

    [~,iA] = min(MAEs);
    [~,iE] = min(RMSEs);

    fprintf('\n=== RANKING (lower is better) ===\n');
    fprintf('MAE (mg/dL)    : %s is best (%.4f)\n', names{iA}, MAEs(iA));
    fprintf('RMSE (mg/dL)   : %s is best (%.4f)\n', names{iE}, RMSEs(iE));


    fprintf('\nCalibration-in-the-large (GMR = yhat/y ratio, ideal=1.00):\n');
    fprintf('  Paper    : %.3f (95%% CI %.3f–%.3f)\n', Mp.GMR,  Mp.GMR_lo95,  Mp.GMR_hi95);
    fprintf('  PaperCal : %.3f (95%% CI %.3f–%.3f)\n', Mpc.GMR, Mpc.GMR_lo95, Mpc.GMR_hi95);
    fprintf('  Ours     : %.3f (95%% CI %.3f–%.3f)\n', Mo.GMR,  Mo.GMR_lo95,  Mo.GMR_hi95);
end

function print_yband_mape(y_obs, y_paperCal, y_ours, idxCommon)
    yb = y_obs(idxCommon);
    edgesY = [0.2 0.5 0.7 1.0 1.5 inf];
    Labs = ["≤0.5","0.5–0.7","0.7–1.0","1.0–1.5",">1.5"];

    binY_sub = discretize(yb, edgesY);
    binY_all = nan(size(y_obs));
    binY_all(idxCommon) = binY_sub;

    fprintf('\nMAPE by observed-y bands:\n');
    for k = 1:numel(Labs)
        mk = idxCommon & (binY_all == k);
        if any(mk)
            m_pCal = mean(abs((y_paperCal(mk) - y_obs(mk))./y_obs(mk))) * 100;
            m_ours = mean(abs((y_ours(mk)     - y_obs(mk))./y_obs(mk))) * 100;
            fprintf('  Y-band %-7s | N=%3d | MAPE paperCal=%.1f%% | ours=%.1f%%\n', ...
                Labs(k), sum(mk), m_pCal, m_ours);
        end
    end
end

function [scale_hat, ci_sc] = global_scale_via_logratio(y_obs, GA, PNA, paper_mu)
% Estimate single multiplicative scale between PAPER and OBS, with CI.
    m  = isfinite(y_obs) & y_obs>0 & isfinite(GA) & isfinite(PNA);
    yo = double(y_obs(m));
    mu = paper_mu(double(GA(m)), double(PNA(m)));
    mu = max(mu, eps);

    r = log( yo ./ mu );
    muR = mean(r); seR = std(r)/sqrt(numel(r));
    ci  = muR + [-1 1]*1.96*seR;

    scale_hat = exp(muR);
    ci_sc     = exp(ci);
end

function T = empirical_scale_pnaBands(tab, y_obs, paper_mu, edges, gaWin)
% Per PNA band: medians of OBS and PAPER & their ratio (+ bootstrap CI).
% Also returns same restricted to GA in [gaWin(1), gaWin(2)]
    GA = double(tab.GA); PNA = double(tab.PNA);
    K  = numel(edges)-1;

    rng(1); B = 1000;

    mkrow = @(lab, idx) band_row(lab, idx, GA, PNA, y_obs, paper_mu, B);

    Rows_all = cell(K,1); Rows_ga = cell(K,1);
    for k=1:K
        inBand = PNA>=edges(k) & PNA<edges(k+1) & isfinite(GA) & isfinite(PNA) ...
               & isfinite(y_obs) & y_obs>0;
        Rows_all{k} = mkrow(sprintf('PNA %g–%g', edges(k), edges(k+1)), inBand);

        inBandGA = inBand & GA>=gaWin(1) & GA<=gaWin(2);
        Rows_ga{k} = mkrow(sprintf('PNA %g–%g | GA %g–%g', edges(k), edges(k+1), gaWin(1), gaWin(2)), inBandGA);
    end

    T.allGA    = vertcat(Rows_all{:});
    T.gaWindow = vertcat(Rows_ga{:});
end

function rowT = band_row(label, idx, GA, PNA, y_obs, paper_mu, B)
    yo = double(y_obs(idx));
    yP = paper_mu(GA(idx), PNA(idx));
    yo = yo(isfinite(yo)); yP = yP(isfinite(yP));
    n = min(numel(yo), numel(yP));

    if n==0
        rowT = table(string(label), 0, NaN, NaN, NaN, NaN, NaN, ...
            'VariableNames', {'Band','N','Med_obs','Med_paper','Ratio','Ratio_lo95','Ratio_hi95'});
        return;
    end

    yo = yo(1:n); yP = yP(1:n);
    med_obs   = median(yo,'omitnan');
    med_paper = median(yP,'omitnan');
    ratio     = med_obs / max(med_paper, eps);

    % bootstrap CI for ratio of medians
    aBoot = nan(B,1);
    for b=1:B
        ii = randi(n, n, 1);
        aBoot(b) = median(yo(ii),'omitnan') / max(median(yP(ii),'omitnan'),eps);
    end
    ci = quantile(aBoot, [0.025 0.975]);

    rowT = table(string(label), n, med_obs, med_paper, ratio, ci(1), ci(2), ...
        'VariableNames', {'Band','N','Med_obs','Med_paper','Ratio','Ratio_lo95','Ratio_hi95'});
end
function [colname, y_obs] = getObservedSCr(T)
    vars = T.Properties.VariableNames; lv = lower(vars);
    cand = vars( contains(lv,'crea') & ~strcmp(lv,'creatinine') );
    if isempty(cand), error('Observed creatinine column not found (looking for "crea__units_*").'); end
    score = zeros(numel(cand),1);
    for i=1:numel(cand)
        li = lower(cand{i});
        score(i) = 10*contains(li,'units') + 3*contains(li,'mg') + 2*contains(li,'umol');
    end
    [~,ix] = max(score); colname = cand{ix};
    y_obs = double(T.(colname));
    li = lower(colname);
    % µmol/L → mg/dL if name suggests µmol
    if contains(li,'umol') || contains(li,'µmol') || contains(li,'umol_l') || contains(li,'umoll')
        y_obs = y_obs / 88.4;
    end
    % guard: drop nonpositive or NaN (will be handled by masks later too)
    bad = ~isfinite(y_obs) | y_obs<=0;
    if any(bad), y_obs(bad) = NaN; end
end


function yhat = predict_ours(tab, flag, models)
% Population-level predictions from your LME (no subject random effects).
    n = height(tab); yhat = nan(n,1); mdl = models{flag + 1};
    for i=1:n
        row = tab(i,:);
        if ~ismember('BW',row.Properties.VariableNames) || isempty(row.BW),        row.BW = 2500 + 700*randn(1,1); end
        if ~ismember('Gender',row.Properties.VariableNames) || isempty(row.Gender),row.Gender = "M"; end
        if ~ismember('cooling',row.Properties.VariableNames) || isempty(row.cooling),row.cooling = 0; end
        if ~ismember('aki',row.Properties.VariableNames) || isempty(row.aki),      row.aki = 0; end
        if ~ismember('Ht',row.Properties.VariableNames) || isempty(row.Ht),        row.Ht = 45 + 4*randn(1,1); end
        if ~ismember('PtID',row.Properties.VariableNames) || isempty(row.PtID),    row.PtID = string(i); end
        Xrow = prepare_data_baby(row, flag);
        logCr = predict(mdl, Xrow, 'Conditional', false); % fixed effects only

        Cr_tot= exp(logCr);

        Cr = min(max(Cr_tot, 0), 4.5);
        yhat(i) = Cr;
    end
end
