%% Neonatal Amikacin PK Simulation and Dosing Optimization
% Chiara Romano & Michela Prunella

clearvars; clc; close all;
rng(0);

%% 1) Virtual Population Generation

simResults=createVirtualPop();

%% 2) Report on VP credibility

T = simResults.neonates;

% Quick affinity checks 
co_affinity_quick_stats(T);

% Hsu Table 2 comparison
hsu = hsu_table2_compare(T); 

% Walther-like stratification 
walther = walther_bw_strat(T, 'co_walther_like.tex'); 

% PDA schematic panels
plot_pda_effects(T, 'pda_effect_schematic');

% Elstad HR under TH plots
elstad_hr_th_plots(T);


%% 3) Dataset selection (no disease, aki, aki and TH)

[dataset, flag] = selectPipelineFromGUI();

%% 4) Fit Model 

[bestFormula_char, lmeBest, rmse, R2] = fitCreatinineModel(dataset,flag);
modelstr = strcat('LMEModel_flag',string(flag),'.mat');
save_path = './SavedModels';
mkdir(save_path)
save(fullfile(save_path, modelstr), 'lmeBest');

%% Find Creatinine models for different categories of newborns 
% for flag = 0:3
%     dataset = getDataset(flag); % or however you select
%     [~, lmeModel, ~, ~] = fitCreatinineModel(dataset, flag);
%     modelstr = strcat('LMEModel_flag',string(flag),'.mat');
%     save_path = './SavedModels';
%     mkdir(save_path)
%     save(fullfile(save_path, modelstr), 'lmeModel');
% end

%% 4) Load Creatinine models
models = cell(4,1);
for i = 0:3
    tmp = load(sprintf('./SavedModels/LMEModel_flag%d.mat', i));
    models{i+1} = tmp.lmeBest;
end

%% 5) Load SimBiology Model
sbioProjObj = sbioloadproject('PBPK_PD_EcoEvolutive.sbproj');
modelObj    = sbioProjObj.m1;

%% 6) REAL POPULATION: Predict Cr and Compute eGFR (mL/min/kg)
R=size(dataset,1);

rng(42);
EstimateCr = zeros(R,1); eGFR = zeros(R,1); CO = zeros(R,1);
Vd_amikacin = zeros(R,1);
if flag == 0
    dataset.DOSE = str2double(string(dataset.DOSE));
end

for i = 1:R
    row = dataset(i,:);

    fprintf("Row %d/%d\n",i,R)
    BW_kg = row.BW/1000;
    if ~isfinite(BW_kg) || BW_kg<=0
            if ismember('BW',row.Properties.VariableNames) && isfinite(row.BW(1))
                BW = row.BW(1);
            elseif ismember('CW',row.Properties.VariableNames) && isfinite(row.CW(1))
                BW = row.CW(1);
            else
                BW = 3.63;
            end
    end
    row.BW = BW;

    if flag == 0 || flag == 1
        [Ht, Vd]  = deriveHeight_VolumeDistribution(row.PNA,row.GA,row.BW);
    elseif flag == 2 || flag ==3
        [Ht, Vd]   = deriveHeight_VolumeDistribution(row.PNA,row.GA,row.BW,row.Gender);
    end

    if flag == 0    % AMICREA-CI
        predict_from = [row.PNA;row.GA;row.BW;row.Ino;row.Ibu];
    else
        predict_from = [row.PNA;row.GA;row.cooling];
    end

    model = RNN_models{flag + 1};


    GA  = row.GA;   % weeks
    pna = row.PNA;  % days

    if GA <= 28                    % like PMA ≤29 when PNA is small
        if      pna < 7,   Vd_per_kg = 0.414;   % their subgroup 1 mean
        elseif  pna <= 28, Vd_per_kg = 0.472;   % subgroup 2 mean
        else               Vd_per_kg = 0.353;   % subgroup 3 mean (tiny n)
        end

    elseif GA <= 33                % like PMA 30–34
        if      pna < 7,   Vd_per_kg = 0.462;   % subgroup 4 mean
        else               Vd_per_kg = 0.454;   % subgroup 5 mean
        end

    else                           % ga ≥ 34  (-> PMA ≥35, all PNA in paper)
        Vd_per_kg = 0.433;                    % subgroup 6 mean
    end


    BW_kg = row.BW/1000;
    Vd_amikacin(i) = Vd_per_kg .* BW_kg;
    row.Ht=Ht;


    Cr = predict(model, predict_from);
    Cr = min(max(Cr, 0), 4.5);
    EstimateCr(i) = Cr;

    k = (GA<37)*0.33 + (GA>=37)*0.45;
    eGFR_std     = k * Ht / Cr;              % mL/min/1.73m2


    % Sistemo BSA con Haycock per coerenza con il resto del codice e
    % commento Karel.
    BSA = 0.024265 .* (BW_kg .^ 0.5378) .* (Ht .^ 0.3964);  % Haycock formula (m^2)


    CL_tot_mL_min= eGFR_std * (BSA/1.73);
    eGFR(i)      = CL_tot_mL_min / BW_kg;    % mL/min/kg

    % HR, SV e CO calcolati ancora nella vecchia maniera. 
    % HR           = 140 + 10.*randn(1);      
    % HR = min(max(HR,100),180); %legato ai battiti medi del bambino
    % SV_mL_per_kg = 2.2 + 0.3.*randn(1);      
    % SV_mL_per_kg = min(max(SV_mL_per_kg,1),4);
    % CO_L_per_min = (HR .* (SV_mL_per_kg .* BW_kg)) / 1000;
    % CO(i) = CO_L_per_min;

    % Cacolo corretto di HR, SV, CO (con altre info tipo TH, PDA, etc..)
    hemi = computeHemodynamicParamsRow(row, flag, 'DerivePDA', true, 'Seed', []);  % [] keeps global RNG
    dataset.HR(i)            = hemi.HR;
    dataset.SV_mL(i)         = hemi.SV_mL;
    dataset.CO_L_min(i) = hemi.CO_L_min;
    dataset.CI_L_min_m2(i)   = hemi.CI_L_min_m2;
    % Optional:
    dataset.PDA(i)           = hemi.PDA;
    dataset.PDA_sev(i)       = hemi.PDA_sev;
    
   
end
dataset.Creatinine  = EstimateCr;
dataset.eGFR = eGFR;
dataset.Vd_amikacin = Vd_amikacin;
writetable(dataset, fullfile(outdir, 'AMICREA_TH_withPhysiology.xlsx'));
save(fullfile(outdir, 'AMICREA_TH_withPhysiology.mat'),'dataset');

%% 7) REAL POPULATION: AFE / AAFE summary  (MECHANISTIC, top-down validation of the Cr→eGFR→CL→Cb pipeline)
% Questo risponde alla nostra domanda: If I trust my Cr→eGFR and physiology
% rules, do the PBPK predictions match observed amikacin levels?" Si! AFE
% vicino ad uno

outdir = './Fitting_Amikacin_prediction';
if ~exist(outdir,'dir'), mkdir(outdir); end

if iscell(dataset.DOSE)
    dataset.DOSE = str2double(string(dataset.DOSE));
end
if iscell(dataset.TIME)
    dataset.TIME = str2double(string(dataset.TIME));
end
if iscell(dataset.CONC)
    dataset.CONC = str2double(string(dataset.CONC));
end

hasPK   = ismember('CONC', dataset.Properties.VariableNames) ...
          && any(isfinite(dataset.CONC));
hasDose = ismember('DOSE', dataset.Properties.VariableNames) ...
          && any(isfinite(dataset.DOSE));
hasTime = ismember('TIME', dataset.Properties.VariableNames) ...
          && any(isfinite(dataset.TIME));

if ~(hasPK && hasDose && hasTime)
    warning('PK evaluation skipped: need CONC, DOSE and TIME with finite values.');
else
    % ---------- (A) TV PBPK predictions driven by eGFR-from-Cr -----------
    fprintf('\n[PK] PBPK-from-Cr: evaluating TV predictions vs observed CONC ...\n');

    fprintf('\n[PK] PBPK-from-Cr: evaluating TIME-VARYING predictions...\n');

    sbioProjObj = sbioloadproject('PBPK_PD_EcoEvolutive.sbproj');
    modelObj    = sbioProjObj.m1;

    report_tv     = evaluateAmikacinAccuracy(dataset, modelObj, true);

    
    % --- per-subject (table) ---
    fprintf('\n[PK] PBPK-from-Cr: evaluating STATIC predictions...\n');
     sbioProjObj = sbioloadproject('PBPK_PD_EcoEvolutive.sbproj');
    modelObj    = sbioProjObj.m1;
    report_static = evaluateAmikacinAccuracy(dataset, modelObj, false);
    
    
    % Copies so we can rename columns
    rs = report_static; rt = report_tv;
    
    % Make sure the key has the same type
    rs.PtID = string(rs.PtID);
    rt.PtID = string(rt.PtID);
    
    % Add suffixes manually
    v = rs.Properties.VariableNames;
    v(2:end) = strcat(v(2:end), '_static');
    rs.Properties.VariableNames = v;
    
    v = rt.Properties.VariableNames;
    v(2:end) = strcat(v(2:end), '_tv');
    rt.Properties.VariableNames = v;
    
    % Full outer join (version-compatible)
    cmp = outerjoin(rs, rt, 'Keys','PtID', 'MergeKeys', true, 'Type', 'full');
    
        
    % Optional: how many obs pairs used per PtID
    pid = string(dataset.PtID);
    good = isfinite(dataset.CONC) & dataset.CONC>0 & isfinite(dataset.TIME);
    pid_g = pid(good);
    [uid, ~, ic] = unique(report_tv.PtID,'stable');
    counts = accumarray( ic( ismember(uid, report_tv.PtID) ), 0); 
    
    writetable(cmp, fullfile(outdir, 'PK_AFE_AAFE_byPt_static_vs_tv.csv'));
    disp(cmp(1:min(10,height(cmp)), :))  % quick peek
    
    % --- pooled (headline numbers) ---
    pooled_static = computePooledFoldMetrics(dataset, modelObj, false);
    pooled_tv     = computePooledFoldMetrics(dataset, modelObj, true);
    
    fprintf('\nPOOLED (all obs–pred pairs)\n');
    fprintf('  STATIC: AAFE=%.3f  AFE=%.3f  P2=%.1f%%  (N=%d)\n', ...
        pooled_static.AAFE, pooled_static.AFE, pooled_static.P2, pooled_static.Npairs);
    fprintf('    TV  : AAFE=%.3f  AFE=%.3f  P2=%.1f%%  (N=%d)\n', ...
        pooled_tv.AAFE, pooled_tv.AFE, pooled_tv.P2, pooled_tv.Npairs);

end

%% 8) Model comparision between our Creatinine models and that in Krzyzanski et al. 2021

run_comparison_paper(models);

%% === PD analysis
% 1) Real dataset with sCr to read complete therapy history
% 2) Read calibration matrix for PD parameters

% Questo risponde alla nostra domanda: Quale è l'effetto PD dello
% scheduling applicato?
outdir = './Impact_Standard_on_PD';
load(".\SavedOptimization\AMICREA_TH_withCalibration.mat") % Read PD patient-specific paramters 
load(fullfile(outdir, 'AMICREA_TH_withPhysiology.mat'))
useTV = true;

if ~exist(outdir,'dir'), mkdir(outdir); end
  
if iscell(dataset.DOSE)
    dataset.DOSE = str2double(string(dataset.DOSE));
end
if iscell(dataset.TIME)
    dataset.TIME = str2double(string(dataset.TIME));
end
if iscell(dataset.CONC)
    dataset.CONC = str2double(string(dataset.CONC));
end

hasPK   = ismember('CONC', dataset.Properties.VariableNames) ...
          && any(isfinite(dataset.CONC));
hasDose = ismember('DOSE', dataset.Properties.VariableNames) ...
          && any(isfinite(dataset.DOSE));
hasTime = ismember('TIME', dataset.Properties.VariableNames) ...
          && any(isfinite(dataset.TIME));

if ~(hasPK && hasDose && hasTime)
    warning('PK evaluation skipped: need CONC, DOSE and TIME with finite values.');
else
    % ---------- (A) TV PBPK predictions driven by eGFR-from-Cr -----------
    fprintf('\n[PK] PBPK-from-Cr: evaluating TV predictions vs observed CONC ...\n');

    fprintf('\n[PK] PBPK-from-Cr: evaluating TIME-VARYING predictions...\n');

    sbioProjObj = sbioloadproject('PBPK_PD_EcoEvolutive.sbproj');
    modelObj    = sbioProjObj.m1;

   
    PtID = dataset.PtID;
    if iscell(PtID), PtID = string(PtID); end
    if iscategorical(PtID), PtID = string(PtID); end
    [IDs, ~, ~] = unique(PtID, 'stable');

    N = numel(IDs);

    doseTargetDefault = "Serum_blood.Aminoglycoside";
    respNameDefault   = "Serum_blood.Cb";
    NrStatus          = "Nr";
    NsStatus          = "Ns";

    dataset.NrSim = zeros(height(dataset),1); dataset.NsSim = zeros(height(dataset),1); dataset.USim = zeros(height(dataset),1);

    for i = 1:N %setdiff(1:N, [212,298,338,387,515,689,698])

        % Read PD parameters
        row   = T_first(i,:);

        % Push physiology + PD theta into model (no dose, no sim)
        theta(1) = row.EC50s0;
        theta(2) = row.EC50r0;
        theta(3) = row.c;
        theta(4) = row.rS;
        theta(5) = row.Emax;
        theta(6) = row.rR;
        theta(7) = row.ku;
        theta(8) = row.sigma;
        theta(9) = row.d;
        p = row.physiology;
    
        apply_calibration_params_only(modelObj, row, p, theta);
        try 
            delete(getdose(modelObj)); 
        end 

        id  = IDs(i);
        sub = dataset(PtID==id, :);
        sub.NrSim = zeros(height(sub),1); sub.NsSim = zeros(height(sub),1); sub.USim = zeros(height(sub),1); 

        % ---- BW (kg) ----
        BW_kg = firstFiniteOr(sub,'BW_kg_',NaN);
        if ~isfinite(BW_kg) || BW_kg<=0
            if ismember('BW',sub.Properties.VariableNames) && isfinite(sub.BW(1))
                BW_kg = sub.BW(1)/1000;
            elseif ismember('CW',sub.Properties.VariableNames) && isfinite(sub.CW(1))
                BW_kg = sub.CW(1)/1000;
            else
                BW_kg = 3.63;
            end
        end

        % ---- eGFR (mL/min/kg) and CO (L/min) ----
        eGFRkg   = firstFiniteOr(sub,'eGFR',NaN);
        CO_L_min = firstFiniteOr(sub,'CO_L_min',NaN);
        if ~isfinite(eGFRkg) || eGFRkg<=0 || ~isfinite(CO_L_min) || CO_L_min<=0
            fprintf('Pt %s: skip (eGFR/CO invalid)\n', string(id));
            continue;
        end

        % ---- GA/PNA for Vd rules ----
        GA  = firstFiniteOr(sub,'GA',NaN);
        PNA = firstFiniteOr(sub,'PNA',NaN);

        % ---- Target Vd (L) ----
        Vd_per_kg = pickVdPerKg(GA,PNA);
        Vd_target = Vd_per_kg * BW_kg;

        % ---- Capacities (L) ----
        BV_per_kg = (GA<37)*0.106 + (GA>=37)*0.080;
        Vc       = BW_kg * BV_per_kg;
        VHearth  = 0.007 * BW_kg;
        Vkidney  = 0.010 * BW_kg;
        Vrem     = max(BW_kg - VHearth - Vkidney, 0);

        % ---- Baseline flows (L/h) ----
        CO_L_h = CO_L_min * 60;
        fh=0.05; fk=0.20; fr=1-fh-fk;
        Qh = max(CO_L_h*fh,1e-6);
        Qk = max(CO_L_h*fk,1e-6);
        Qr = max(CO_L_h*fr,1e-6);
        s  = CO_L_h/(Qh+Qk+Qr);  Qh = Qh*s; Qk = Qk*s; Qr = Qr*s;

        % ---- Baseline CL (L/h) ----
        CL_L_h = eGFRkg * BW_kg * 0.06;

        % ---- Dose & Obs ----
        dAmt = sub.DOSE;  if ~isnumeric(dAmt), dAmt = str2double(string(dAmt)); end
        dTim = sub.TIME;
        mDose= isfinite(dAmt) & isfinite(dTim);
        dAmt = dAmt(mDose); dTim = dTim(mDose);

        mObs = isfinite(sub.CONC) & isfinite(sub.TIME);
        tObs = sub.TIME(mObs);  yObs = sub.CONC(mObs);
        if isempty(tObs)
            fprintf('Pt %s: skip (no observations)\n', string(id));
            continue;
        end
        tStop = max(tObs) + 24;

        % ===== model clone =====
        m = tryCopyModel(modelObj);
        cs = getconfigset(m);
        cs.StopTime = tStop;

        apply_calibration_params_only(m, row, p, theta);
        % ---- Capacities ----
        setCapL(m,'Serum_blood', Vc);
        setCapL(m,'Hearth',      VHearth);   % tolerate "Hearth"
        setCapL(m,'Heart',       VHearth);   % and "Heart"
        setCapL(m,'Kidney',      Vkidney);
        setCapL(m,'Remainder',   Vrem);

        % ---- Match Vd via Kp remainder ----
        Kp_h = getParamOr(m,'PartitionCoefficientHearth',  0.25);
        if isnan(Kp_h), Kp_h = getParamOr(m,'PartitionCoefficientHeart', 0.25); end
        Kp_k = getParamOr(m,'PartitionCoefficientKidney',  0.30);
        num  = Vd_target - Vc - Kp_h*VHearth - Kp_k*Vkidney;
        Kp_r = max(num / max(Vrem,eps), 0.05);
        setParamIfExists(m,'PartitionCoefficientHearth',     Kp_h);
        setParamIfExists(m,'PartitionCoefficientHeart',      Kp_h);
        setParamIfExists(m,'PartitionCoefficientKidney',     Kp_k);
        setParamIfExists(m,'PartitionCoefficientRemainder',  Kp_r);

        % ---- Push baseline flows/CL (DO NOT also set EliminationRate unless needed) ----
        setParamIfExists(m,'BloodFlowRateHearth',     Qh);
        setParamIfExists(m,'BloodFlowRateHeart',      Qh);
        setParamIfExists(m,'BloodFlowRateKidney',     Qk);
        setParamIfExists(m,'BloodFlowRateRemainder',  Qr);
        setParamIfExists(m,'CL',                      CL_L_h);
        % setParamIfExists(m,'EliminationRate', CL_L_h/max(Vc,1e-9));  % <- leave OFF to avoid double-clearance
        setParamIfExists(m,'TBW', BW_kg);
        setParamIfExists(m,'GFR', eGFRkg);

        % ---- Dose target ----
        doseTarget = resolveDoseTarget(m, doseTargetDefault);
        if doseTarget == ""
            fprintf('Pt %s: skip (no dose target)\n', string(id));
            continue;
        end

        % ---- TV physiology events ----
        tch = [];
        if useTV
            [nEvents, tch] = addTimeVaryingPhysiology(m, sub, BW_kg, Vc);
            fprintf('Pt %s: events added = %d\n', string(id), nEvents);
        end

        % ---- Make solver hit obs & event times ----
        if ~isempty(tObs)
            cs.SolverOptions.OutputTimes = unique([0; double(tObs(:)); double(tch(:)); cs.StopTime]);
        else
            cs.SolverOptions.OutputTimes = unique([0; double(tch(:)); cs.StopTime]);
        end
        cs.RuntimeOptions.StatesToLog = {'Serum_blood.Cb','Serum_blood.Ns','Serum_blood.Nr','u'};

        % ---- Simulate ----
        try
            if ~isempty(dTim)
                d = sbiodose('Amikacin','schedule');
                d.TargetName  = char(doseTarget);
                d.Amount      = dAmt(:);
                d.Time        = dTim(:);
                d.AmountUnits = 'milligram';
                d.Rate        = d.Amount ./ 1;
                d.RateUnits   = 'milligram/hour';
                d.TimeUnits   = 'hour';
                d.Active      = true;

                m.Species(3,1).InitialAmount = 8.8087e+05;
                m.Species(4,1).InitialAmount = 2.2194e+05;
                sd = sbiosimulate(m, d);
            else
                m.Species(3,1).InitialAmount = 8.8087e+05;
                m.Species(4,1).InitialAmount = 2.2194e+05;
                sd = sbiosimulate(m);
            end
        catch ME
            fprintf('Pt %s: simulation failed: %s\n', string(id), ME.message);
            continue;
        end

        % ---- Response extraction ----
        dataNames = string(sd.DataNames);
        idx = find(dataNames == respNameDefault, 1);
        idxnr = find(dataNames == NrStatus, 1);
        idxns = find(dataNames == NsStatus, 1);
        idxu = find(dataNames == 'u', 1);
        if isempty(idx), idx = find(dataNames == "Cb", 1); end
        if isempty(idx)
            fprintf('Pt %s: response not found\n', string(id));
            continue;
        end
        
        tSim = double(sd.Time(:));

        for t=1:height(sub)
            sub.NrSim(t) =  safeInterp(tSim, double(sd.Data(:, idxnr)), double(sub.TIME(t)));
            sub.NsSim(t) =  safeInterp(tSim, double(sd.Data(:, idxns)), double(sub.TIME(t)));
            sub.USim(t) =  safeInterp(tSim, double(sd.Data(:, idxu)), double(sub.TIME(t)));
        end

        dataset(PtID==id, :) = sub;
    end
    writetable(dataset, fullfile(outdir, 'PD_DynamicsStandardTherapy_AMICREA_TH.xlsx'));
end

%% Plot Residuals (questa parte non l'ho verificata)
% timeVec = data.(data.Properties.IndependentVariableName);
% obsData = data.CONC;
% indicesToKeep = ~isnan(obsData);
% timeVec = timeVec(indicesToKeep);
% 
% % Get the residuals from the fitting results.
% indRes = nlmeResults.stats.ires(indicesToKeep);
% popRes = nlmeResults.stats.pres(indicesToKeep);
% 
% % Plot the residuals. Get a handle to the plot to be able to add more data
% % to it later.
% resplot = figure;
% plot(timeVec,indRes,'d','MarkerFaceColor','b','markerEdgeColor','b');
% hold on;
% plot(timeVec,popRes,'d','MarkerFaceColor','w','markerEdgeColor','b');
% hold off;
% 
% % Create a reference line representing a zero residual, and set its
% % properties to omit this line from the plot legend.
% refl = refline(0,0);
% refl.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% % Label the plot.
% title('Residuals versus Time');
% xlabel('Time');
% ylabel('Residuals');
% legend({'Individual','Population'});

%% 9) VIRTUAL POPULATION: Predict Cr and Compute eGFR (mL/min/kg) 
N=size(simResults.neonates,1);
Estimate = zeros(N,1);
eGFR     = zeros(N,1);
CO = zeros(N,1);
Cr_tot = zeros(N,1);

for i = 1:N
    row = simResults.neonates(i,:);

    hasAKI = row.aki == 1;
    isCooled = row.cooling == 1;

    % Determine flag
    if ~hasAKI && ~isCooled
        local_flag = 0;
    elseif hasAKI && ~isCooled
        local_flag = 1;
    elseif ~hasAKI && isCooled
        local_flag = 2;
    else
        local_flag = 1;
    end

    model = models{local_flag + 1};


    GA  = row.GA;
    BW_kg = row.BW/1000;
    Ht   = row.Ht;

    logCr = predict(model, prepare_data_baby(row,local_flag));
    Cr_tot(i)= exp(logCr);
    Cr = min(max(Cr_tot(i), 0), 4.5);

    Estimate(i) = Cr;

    k = (GA<37)*0.33 + (GA>=37)*0.45;
    eGFR_std     = k * Ht / Cr;              % mL/min/1.73m2
    BSA          = row.BSA;  
    CL_tot_mL_min= eGFR_std * (BSA/1.73);
    eGFR(i)      = CL_tot_mL_min / BW_kg;    % mL/min/kg

    
    CO_L_per_min = row.CO_L_min;
    CO(i)=CO_L_per_min;

end
simResults.neonates.Creatinine  = Estimate;
simResults.neonates.eGFR = eGFR;
simResults.neonates.CardiacOutput = CO;

% ---- Remove rows with saturated creatinine (== 4.5 mg/dL)
mask = simResults.neonates.Creatinine < 4.5;
simResults.neonates = simResults.neonates(mask,:);


% ---- Remove rows with saturated creatinine (== 2.5 mg/dL)
% mask2 = simResults.neonates.Creatinine > 0.2;
% simResults.neonates = simResults.neonates(mask2,:);

fprintf('Removed %d rows with Cr >= 4.5 mg/dL (saturating value).\n', sum(~mask));
%fprintf('Removed %d rows with Cr < 0.2 mg/dL (saturating value).\n', sum(~mask2));
mask = simResults.neonates.Creatinine > 0.2;
simResults.neonates = simResults.neonates(mask,:);

fprintf('Removed %d rows with Cr < 0.2 mg/dL (saturating value).\n', sum(~mask));

%% 10) PLOT to check Creatinine consistency
% Scatter plot of all individual data points
figure
scatter(simResults.neonates.PNA, simResults.neonates.Creatinine, 'b.'); 
hold on;

% Compute mean Creatinine per PNA
[PNA_unique, ~, idx] = unique(simResults.neonates.PNA);
meanCr = accumarray(idx, simResults.neonates.Creatinine, [], @mean);

% Plot the mean trend line
plot(PNA_unique, meanCr, 'r-', 'LineWidth', 2);

xlabel('PNA (days)');
ylabel('Creatinine (mg/dL)');
title('Creatinine vs. PNA with Mean Trend');
legend('Data points', 'Mean Creatinine');
grid on;

idxOutliers = find(Cr_tot > 2.5);
outlierRows = simResults.neonates(idxOutliers, :);
disp(Cr_tot(idxOutliers))
disp(outlierRows);

%% 11) PLOT to check Creatinine consistency (pannello inferiore della figura 4 report)
% Scatter plot of all individual data points
figure
scatter(simResults.neonates.PNA, simResults.neonates.Creatinine, 'b.'); 
hold on;

% Compute centiles as you already do ...
PNA_unique = unique(simResults.neonates.PNA);
PNA_unique = PNA_unique(~isnan(PNA_unique));

pgrid = [0.10 0.25 0.50 0.75 0.90 0.95];
cent  = nan(numel(PNA_unique), numel(pgrid));
n_per_day = zeros(numel(PNA_unique),1);

for i = 1:numel(PNA_unique)
    vals = simResults.neonates.Creatinine(simResults.neonates.PNA == PNA_unique(i));
    vals = vals(~isnan(vals));
    n_per_day(i) = numel(vals);
    if ~isempty(vals)
        cent(i,:) = quantile(vals, pgrid);
    end
end

% Optional: mask days with very small N (e.g., < 5) to avoid degenerate quantiles
mask = n_per_day >= 5;
x = PNA_unique(mask);
c = cent(mask,:);

figure('Color','w','Position',[100 100 700 420]); hold on;
scatter(simResults.neonates.PNA, simResults.neonates.Creatinine, 10, 'MarkerEdgeAlpha', 0.27, 'DisplayName','Data');

% plot each centile separately so legend shows both
plot(x, c(:,3), '-',  'LineWidth', 2,   'DisplayName','p50');
plot(x, c(:,2), '--', 'LineWidth', 1.5, 'DisplayName','p25');
plot(x, c(:,4), '--', 'LineWidth', 1.5, 'DisplayName','p75');
plot(x, c(:,1), ':',  'LineWidth', 1.2, 'DisplayName','p10');
plot(x, c(:,5), ':',  'LineWidth', 1.2, 'DisplayName','p90');
plot(x, c(:,6), '-.', 'LineWidth', 1.2, 'DisplayName','p95');

xlabel('PNA (days)'); ylabel('Serum creatinine (mg/dL)');
title('Virtual population: Serum creatinine centiles vs PNA');
legend('Location','best');
grid on; xlim([0 max(PNA_unique)]);
ylim([0 2])
xlim([0 30])
box on;

% Quick diagnostic: are p25 and p75 identical anywhere?
d = c(:,4) - c(:,2);
fprintf('p75 - p25: min=%.4g, median=%.4g, max=%.4g (over days with N>=5)\n', ...
        min(d,[],'omitnan'), median(d,'omitnan'), max(d,[],'omitnan'));


%% Effective calibration with Nielsen model after amikacin adjustments wrt genta
% row = simResults.neonates(1,:);
% 
% baseDir  = '.\Bacterial_model_calibration';   % where g1_25.csv etc live
% schedule = 'g2_25';                           % <- pick one (GA<32 in your naming)
% optsCal  = struct('duration_h', 72, 'fu', 0.9, 'kf', 1, 'equiv_factor', 4);
% 
% [theta_hat, sim, p, meta] = pbpk_calibrate_from_genta(row, schedule, baseDir, optsCal);
% 
% fprintf('Equivalence used: gent %.1f mg/kg q%d h  ->  amikacin %.1f mg/kg q%d h\n', ...
%     meta.gent_dose_mgkg, meta.tau_h, meta.ami_dose_mgkg, meta.tau_h);
% 
% % quick plot vs target
% targets = load_targets_for_schedule(schedule, baseDir);
% figure('Color','w'); hold on;
% plot(sim.t, log10(sim.y(:,7)+sim.y(:,8)), 'LineWidth', 2);
% scatter(targets.t, log10(targets.N), 35, 'filled');
% xlabel('Time (h)'); 
% ylabel('log_{10} CFU/mL');   % instead of 'Total CFU/mL'
% grid on; box on;
% legend('Amikacin (equiv) model','Gent targets');
% title(sprintf('PD calibration for %s (newborn %s)', schedule, string(row.PtID)));


%% Effective calibration with Nielsen model after amikacin adjustments wrt genta for ALL schedules 
% idx = find(simResults.neonates.GA == 25, 1, 'first');
% row = simResults.neonates(idx, :);
% 
% schedules = {'g1_25','g2_25','g3_25'};  % or _40, etc., matching your baby
% [theta_hat1,p1, k_eq1] = pbpk_calibrate_joint(row, schedules, baseDir, struct('duration_h',72));
% 
% validate_across_schedules(row, baseDir, theta_hat1, k_eq1, schedules, 72); %per check con fig. 4 Nielsen
% 
% 
% idx = find(simResults.neonates.GA == 40, 1, 'first');
% 
% row = simResults.neonates(idx, :);
% 
% schedules = {'g1_40','g2_40','g3_40'};  % or _40, etc., matching your baby
% [theta_hat2,p2, k_eq2] = pbpk_calibrate_joint(row, schedules, baseDir, struct('duration_h',72));
% 
% validate_across_schedules(row, baseDir, theta_hat2, k_eq2, schedules, 72);

%% CHECK CALIBRATION 
warning on; 
sbioProjObj = sbioloadproject('PBPK_PD_EcoEvolutive.sbproj');
modelObj    = sbioProjObj.m1;
load(".\SavedOptimization\DigitalTwins_pop.mat"); % -> will write into simResults the new population (15% of patients with PATH have AKI)
 
N = height(simResults.neonates);
simResults.neonates.EC50s0 = nan(N,1); simResults.neonates.EC50r0 = nan(N,1); simResults.neonates.c = nan(N,1);
simResults.neonates.rS = nan(N,1); simResults.neonates.Emax = nan(N,1); simResults.neonates.rR = nan(N,1);
simResults.neonates.ku = nan(N,1); simResults.neonates.sigma = nan(N,1); simResults.neonates.d = nan(N,1);
simResults.neonates.physiology = repmat(struct('Serum_blood', [], 'Hearth', [], 'Kidney', [],'Remainder', [],'BloodFlowRateHearth', [], 'BloodFlowRateKidney', [],...
    'BloodFlowRateRemainder', [], 'PartitionCoefficientHearth', [], 'PartitionCoefficientKidney', [], 'PartitionCoefficientRemainder', [],'CL', [], ...
    'kf', [], 'fu', [], 'aRS', [], 'aSR', [],'K', []),N,1);
 
for i = 1:N
    row   = simResults.neonates(i,:);
 
    guards = classical_schedules(row);
 
    [theta_hat,p] = pbpk_calibrate_joint(row, struct('duration_h',72, 'guard_schedules',{guards}));
    % EC50s0, EC50r0, c, rS, Emax, rR, ku, sigma, d
    row.EC50s0 = theta_hat(1);
    row.EC50r0 = theta_hat(2);
    row.c = theta_hat(3);
    row.rS = theta_hat(4);
    row.Emax = theta_hat(5);
    row.rR = theta_hat(6);
    row.ku = theta_hat(7);
    row.sigma = theta_hat(8);
    row.d = theta_hat(9);
    row.physiology = p;
    simResults.neonates(i,:) = row;
 
    % plot_classical_schedule_for_baby(modelObj, row, theta_hat, p);
end
save(".\SavedOptimization\DigitalTwins_pop_withCalibration.mat","simResults") 

%% CALIBRATE PD ON REAL POP 
warning on; 
sbioProjObj = sbioloadproject('PBPK_PD_EcoEvolutive.sbproj');
modelObj    = sbioProjObj.m1;

% Considering the first row for each patient
[~, idx] = unique(dataset.PtID, 'first');

% Estrai solo la prima riga per ciascun ID
T_first = dataset(idx, :);
N = height(T_first);
T_first.EC50s0 = nan(N,1); T_first.EC50r0 = nan(N,1); T_first.c = nan(N,1);
T_first.rS = nan(N,1); T_first.Emax = nan(N,1); T_first.rR = nan(N,1);
T_first.ku = nan(N,1); T_first.sigma = nan(N,1); T_first.d = nan(N,1);
T_first.physiology = repmat(struct('Serum_blood', [], 'Hearth', [], 'Kidney', [],'Remainder', [],'BloodFlowRateHearth', [], 'BloodFlowRateKidney', [],...
    'BloodFlowRateRemainder', [], 'PartitionCoefficientHearth', [], 'PartitionCoefficientKidney', [], 'PartitionCoefficientRemainder', [],'CL', [], ...
    'kf', [], 'fu', [], 'aRS', [], 'aSR', [],'K', []),N,1);
 
for i = 1:N
    fprintf("Calibrating PtID: %d\n",i)
    row   = T_first(i,:);
 
    guards = classical_schedules(row);
 
    [theta_hat,p] = pbpk_calibrate_joint(row, struct('duration_h',72, 'guard_schedules',{guards}));
    % EC50s0, EC50r0, c, rS, Emax, rR, ku, sigma, d
    row.EC50s0 = theta_hat(1);
    row.EC50r0 = theta_hat(2);
    row.c = theta_hat(3);
    row.rS = theta_hat(4);
    row.Emax = theta_hat(5);
    row.rR = theta_hat(6);
    row.ku = theta_hat(7);
    row.sigma = theta_hat(8);
    row.d = theta_hat(9);
    row.physiology = p;
    T_first(i,:) = row;
 
    % plot_classical_schedule_for_baby(modelObj, row, theta_hat, p);
end
save(fullfile(outdir, "AMICREA_TH_withCalibration.mat"),"T_first") 

%%  Questa sezione è servita per fare il merge tra la figura ottenuta dal file NONMEM_dataset_analysis e quella nella precedente sezione
% % Open the two saved figures (do not show)
% f1 = openfig('nonmem.fig','invisible');
% f2 = openfig('virtual.fig','invisible');
% 
% % Get the FIRST data axes from each figure (exclude colorbars/uipanels)
% srcAx1 = findall(f1,'Type','axes','-not','Tag','Colorbar');
% srcAx2 = findall(f2,'Type','axes','-not','Tag','Colorbar');
% srcAx1 = srcAx1(1);   % pick primary axes if multiple
% srcAx2 = srcAx2(1);
% 
% % New stacked figure (2 rows, 1 col)
% F = figure('Color','w','Position',[100 100 700 900]);
% tl = tiledlayout(F,2,1,'TileSpacing','compact','Padding','compact');
% 
% % ---- Top tile: NONMEM
% axT = nexttile(tl); hold(axT,'on');
% copyobj(allchild(srcAx1), axT);
% xlim(axT, xlim(srcAx1)); ylim(axT, ylim(srcAx1));
% xlabel(axT, srcAx1.XLabel.String);
% ylabel(axT, srcAx1.YLabel.String);
% title(axT,  get(get(srcAx1,'Title'),'String'));
% grid(axT,'on'); box(axT,'on');
% % Copy legend if present
% lg1 = findobj(f1,'Type','Legend');
% if ~isempty(lg1)
%     legend(axT, lg1.String, 'Location', lg1.Location);
% end
% 
% % ---- Bottom tile: Virtual
% axB = nexttile(tl); hold(axB,'on');
% copyobj(allchild(srcAx2), axB);
% xlim(axB, xlim(srcAx2)); ylim(axB, ylim(srcAx2));
% xlabel(axB, srcAx2.XLabel.String);
% ylabel(axB, srcAx2.YLabel.String);
% title(axB,  get(get(srcAx2,'Title'),'String'));
% grid(axB,'on'); box(axB,'on');
% % Copy legend if present
% lg2 = findobj(f2,'Type','Legend');
% if ~isempty(lg2)
%     legend(axB, lg2.String, 'Location', lg2.Location);
% end
% 
% % Optional: same Y-limits
% yl = [min([ylim(axT) ylim(axB)]) max([ylim(axT) ylim(axB)])];
% ylim(axT,yl); ylim(axB,yl);
% 
% % Clean up source figures
% close(f1); close(f2);
% 

%% Prepare REAL DATASET for OPTIMIZATION
% Working on all columns
N=size(T_first,1); 

load ./SavedOptimization/AMICREA_CI_withCalibration.mat

if flag==0                      % AMICREA-CI
    % AKI, COOLING is absent
    aki=zeros(N,1); 
    cooling=zeros(N,1);
    T_first.cooling = cooling;

    % Gender is absent
    TBW = T_first.BW.*0.755/1000;
    T_first.TBW = TBW;
    T_first.Gender = repmat("m",N,1);
    [Ht, Vd]  = deriveHeight_VolumeDistribution(T_first.PNA,T_first.GA,T_first.BW);

else                            % AMICREA-TH
    % 15% have AKI
    aki=zeros(N,1); 
    idx = randperm(N, round(0.15*N)); % seleziona casualmente il 15% degli indici
    aki(idx) = 1;
    
    % Gender is specified
    TBW = zeros(N,1);
    TBW(strcmp(T_first.Gender, 'm')) = 0.79 * T_first.BW(strcmp(T_first.Gender, 'm'));  % Males: 79% of BW
    TBW(strcmp(T_first.Gender, 'v')) = 0.72 * T_first.BW(strcmp(T_first.Gender, 'v'));  % Females: 72% of BW
    T_first.TBW = TBW/1000;

    [Ht, Vd]  = deriveHeight_VolumeDistribution(T_first.PNA,T_first.GA,T_first.BW,T_first.Gender);
  
end

T_first.aki = aki; 
T_first.Ht = Ht;
T_first.BW_kg_ = T_first.BW./1000;
T_first.BSA = 0.024265 .* (T_first.BW_kg_ .^ 0.5378) .* (T_first.Ht .^ 0.3964); % Haycock formula
    
% Set MIC
micValue  = [2, 3, 4, 6, 8];  % tipici sono 2/4/6. 
micFreq   = [6, 2, 1, 0.5, 0.5];

% Sample MIC values from a discrete distribution using randsample
MIC = randsample(micValue , N, true , micFreq);
T_first.MIC = MIC';
T_first.PtID = string(T_first.PtID);

disp("Adding required fields to real cohort for further optimization")

%% 3) Optimize drug regimen
rng(42);
warning off;
sbioProjObj = sbioloadproject('PBPK_PD_EcoEvolutive.sbproj');
modelObj    = sbioProjObj.m1;

modelObj.Species(3,1).InitialAmount = 8.8087e+05;
modelObj.Species(4,1).InitialAmount = 2.2194e+05;

% load ./Bacterial_model_calibration/DigitalTwins_jointCalibration.mat

load ./SavedModels/RNN_models.mat

% N=size(simResults.neonates,1);
% FOR REAL DATASETS
N=size(T_first,1); 


opt.bounds.DoseAmount_mg_per_kg = [10, 25];

opt.doseName   = 'Serum_blood.Aminoglycoside';
opt.doseUnits  = 'milligram';
opt.isInfusion = true;
opt.infDur     = 1;     % h
opt.rateUnits  = 'milligram/hour';
opt.concName   = 'Cb';
opt.postWindow = 72;    % h

resultsOptPop = table('Size', [0,32], 'VariableNames',{'PtID','GA','PNA','PMA','Gender','BW','Ht','Vd','aki','cooling','Ibu','Ino', ...
                'HR','BSA','SV_mL','CO_L_min','TBW','MIC','Vd_amikacin','Creatinine','eGFR','CardiacOutput','GFR_L_per_min','DoseOpt_mg_kg', ...
                'Dose_mg','Cmax','TCmin','CmaxMIC','AUCMIC','CFUBact','NsLog','NrLog'}, ...
                'VariableTypes',["string",repmat("double",1,3), "string", repmat("double",1,27)]);


target.Cmax_min        = 0;      % non necessario se usi solo il rapporto
target.Cmax_cap_hard   = 60;     % cap duro sicurezza
target.Cmax_cap_soft   = 35;     % penalità soft oltre 35

target.TCmin_max       = 48;     % TCmin deve essere ≤ 48 h

target.CmaxMIC_lb      = 8;      % banda obiettivo del rapporto
target.CmaxMIC_ub      = 10;

% Lasciali ma NON li useremo nei vincoli/obiettivo ora
target.AUC_MIC_lb      = 21.4;
target.AUC_MIC_goal2   = 62.5;

scale_factors=[];
for i = 1:N
   
    resultOpt = table('Size', [0,32], 'VariableNames',{'PtID','GA','PNA','PMA','Gender','BW','Ht','Vd','aki','cooling','Ibu','Ino', ...
                'HR','BSA','SV_mL','CO_L_min','TBW','MIC','Vd_amikacin','Creatinine','eGFR','CardiacOutput','GFR_L_per_min','DoseOpt_mg_kg', ...
                'Dose_mg','Cmax','TCmin','CmaxMIC','AUCMIC','CFUBact','NsLog','NrLog'}, ...
                'VariableTypes',["string",repmat("double",1,3), "string", repmat("double",1,27)]);

    row   = T_first(i,:);

    total_hours_therapy = 0; CFUBact=[0,nan(1,19)]; 
    decline_duration = 0;
    
    % Push physiology + PD theta into model (no dose, no sim)
    theta(1)=row.EC50s0;
    theta(2)=row.EC50r0;
    theta(3)=row.c;
    theta(4)=row.rS;
    theta(5)=row.Emax;
    theta(6)=row.rR;
    theta(7)=row.ku;
    theta(8)=row.sigma;
    theta(9)=row.d;
    p=row.physiology;

    apply_calibration_params_only(modelObj, row, p, theta);
    try 
        delete(getdose(modelObj)); 
    end 

    % --- Prima del while / for loop ---
    tolAbs = 1e-6;   % tolleranza assoluta (mg/L) -> accetta differenze piccolissime
    tolRel = 1e-3;   % tolleranza relativa (0.1%) -> evita falsi allarmi su valori grandi
    
    Cb_end_prev = NaN;   % inizializza

    maxIters  = 20;         % safety cap: ~60 intervals is plenty
    maxHours  = 240;        % hard max total therapy hours

    j = 0;


    while total_hours_therapy <= 144 || decline_duration < 3 
        j=j+1;
        hasAKI   = row.aki == 1;
        isCooled = row.cooling == 1;
    
        if ~hasAKI && ~isCooled    % AMICREA-CI
            local_flag = 0;
            predict_from = [row.PNA;row.GA;row.BW;row.Ino;row.Ibu];
        elseif hasAKI && ~isCooled % CREA-AKITH
            local_flag = 1;
            predict_from = [row.PNA;row.GA;row.aki];
        elseif ~hasAKI && isCooled % AMICREA-TH
            local_flag = 2;
            predict_from = [row.PNA;row.GA;row.cooling];
        else
            local_flag = 1;
            predict_from = [row.PNA;row.GA;row.aki];
        end
    
        model = RNN_models{local_flag + 1};

        GA     = row.GA;
        BW_kg  = row.BW/1000;
        TBW_kg = row.TBW;

        % --- Renal function
        Cr  = predict(model, predict_from);
        Cr_pred = min(max(Cr, 0.4), 4.5); % prima stava a 0.9, VINCOLATA TROPPO secondo me -> il max val of creat in creacool è 4

        Ht     = row.Ht;
        k      = (GA < 37)*0.33 + (GA >= 37)*0.45;
        eGFR_std = k * Ht / Cr_pred;               % mL/min/1.73m^2
        BSA      = row.BSA;        % CORRETTO ora
        CL_tot_mL_min = eGFR_std * (BSA / 1.73);   % mL/min
        eGFRkg  = CL_tot_mL_min / BW_kg;           % mL/min/kg

        % --- Volumes (L)
        BV_per_kg = (GA<37)*0.106 + (GA>=37)*0.080;   % L/kg
        Vc        = BW_kg * BV_per_kg;                % central blood/plasma volume
        Vheart    = 0.007 * BW_kg;                    % ~0.6% BW
        Vkidney   = 0.010 * BW_kg;                    % ~1.0% BW
        Vrem      = max(BW_kg - Vheart - Vkidney, 0);

        % Push capacities (species are concentrations)
        sb=sbioselect(modelObj,'Name','Serum_blood','Type','compartment'); sb.Capacity = Vc;      sb.CapacityUnits='liter';
        ht=sbioselect(modelObj,'Name','Heart','Type','compartment');       ht.Capacity = Vheart;  ht.CapacityUnits='liter';
        kd=sbioselect(modelObj,'Name','Kidney','Type','compartment');      kd.Capacity = Vkidney; kd.CapacityUnits='liter';
        rm=sbioselect(modelObj,'Name','Remainder','Type','compartment');   rm.Capacity = Vrem;    rm.CapacityUnits='liter';

        % --- Whole-body Vd target via Kp (Vss = Vc + Σ Kp*Vt)
        Vd_target = row.Vd_amikacin;   % L
        Kp_heart  = getparamvalue(modelObj,'PartitionCoefficientHeart',  0.25);
        Kp_kidney = getparamvalue(modelObj,'PartitionCoefficientKidney', 0.30);
        num       = Vd_target - Vc - Kp_heart*Vheart - Kp_kidney*Vkidney;
        Kp_rem    = max(num / max(Vrem,eps), 0.05);   % clamp to ≥0.05
        setparam(modelObj,'PartitionCoefficientHeart',      Kp_heart);
        setparam(modelObj,'PartitionCoefficientKidney',     Kp_kidney);
        setparam(modelObj,'PartitionCoefficientRemainder',  Kp_rem);

        % --- Cardiac output -> flows (L/h) that sum to CO
        CO_L_per_h = row.CO_L_min * 60;
        frac_h = 0.05; frac_k = 0.20; frac_r = 1 - frac_h - frac_k;
        Qh = max(CO_L_per_h*frac_h,1e-6);
        Qk = max(CO_L_per_h*frac_k,1e-6);
        Qr = max(CO_L_per_h*frac_r,1e-6);
        % renormalize to sum exactly to CO
        scale = CO_L_per_h/(Qh+Qk+Qr);
        Qh = Qh*scale; Qk = Qk*scale; Qr = Qr*scale;

        setparamLh(modelObj,'BloodFlowRateHeart',     Qh);   % L/h
        setparamLh(modelObj,'BloodFlowRateKidney',    Qk);   % L/h
        setparamLh(modelObj,'BloodFlowRateRemainder', Qr);   % L/h

        % --- Clearance (L/h) and optional kel
        CL_L_per_h = eGFRkg * BW_kg * 0.06;                 % (mL/min/kg) -> L/h
        setparamLh(modelObj,'CL', CL_L_per_h);              % used by: CL * Cb
        kel_per_h = CL_L_per_h / max(Vc,1e-9);              % optional (if referenced)
        setparam(modelObj,'EliminationRate', kel_per_h);

        % TBW & GFR (your original intent)
        tbw=sbioselect(modelObj,'Name','TBW','Type','parameter'); tbw.Value = TBW_kg;
        gfr=sbioselect(modelObj,'Name','GFR','Type','parameter'); gfr.Value = eGFRkg;
        
        % se vogliamo ricalibrare ogni volta basta reinserire la
        % calibrazione qui.. 

        % --- Optimize regimen and simulate
        if j==1
            target.MIC = row.MIC;
            [Xopt,fval] = optimizeRegimen(modelObj, getconfigset(modelObj), target, opt, BW_kg, []);
            [AUC,Cmax,TCmin, IC, simData] = simulateMetrics(modelObj, getconfigset(modelObj), Xopt, opt, BW_kg, []);
            [Cb_start, Cb_end] = get_Cb_start_end(simData);
            Cb_end_prev = Cb_end;
        else
            baseline_MIC = row.MIC;
            u_val = IC(5);
            target.MIC = MIC_from_u(u_val, baseline_MIC);
            [Xopt,fval] = optimizeRegimen(modelObj, getconfigset(modelObj), target, opt, BW_kg, IC);
            [AUC,Cmax,TCmin, IC, simData] = simulateMetrics(modelObj, getconfigset(modelObj), Xopt, opt, BW_kg, IC);
            [Cb_start, Cb_end] = get_Cb_start_end(simData);
            jump = Cb_start - Cb_end_prev;
            if abs(jump) > max(tolAbs, tolRel*max(1e-9,abs(Cb_end_prev)))
                warning('DISCONTINUITA'' Cb: prev_end=%.6g, curr_start=%.6g (jump=%.3g)', Cb_end_prev, Cb_start, jump);
            end
            Cb_end_prev = Cb_end;
        end

        % Extract Ns & Nr (model variable names must match your SimBiology model)
        tend = simData.Time(end);
        Ns_end = simData.Data(end, strcmp(simData.DataNames,'Ns'));
        Nr_end = simData.Data(end, strcmp(simData.DataNames,'Nr'));
        
        % Ensure log10 scale (if your CFUBact is already a log10 metric, keep consistent):
        epsN = 1e-12;
        Ns_log = log10(max(Ns_end, epsN));
        Nr_log = log10(max(Nr_end, epsN));
        
        % Store
        resultOpt.NsLog(j) = Ns_log;
        resultOpt.NrLog(j) = Nr_log;


        % --- stopping logic
        total_hours_therapy = total_hours_therapy + TCmin;
        CFUBact(j) = IC(4);
        if total_hours_therapy > 96
             isNeg = CFUBact(j) < 0;
             decline_duration = isNeg * (decline_duration + 1);   % keeps your consecutive-count
        end

        % --- store dose & outputs
        dose_mg_per_kg = Xopt;
        dose_mg        = dose_mg_per_kg * BW_kg;

        row.GFR_L_per_min   = CL_tot_mL_min/1000;   % L/min
        row.DoseOpt_mg_kg   = dose_mg_per_kg;
        row.Dose_mg         = dose_mg;    
        row.Cmax            = Cmax;
        row.TCmin           = TCmin;
        row.MIC             = target.MIC;
        row.CmaxMIC         = Cmax/target.MIC;
        row.AUCMIC          = AUC/target.MIC;
        row.Creatinine      = Cr_pred;
        row.CFUBact         = CFUBact(j);

        fprintf('[Pt %s]\t[Iter %d]: dose=%.2f mg/kg | MIC=%.2f | Cmax=%.2f | Cmax/MIC=%.2f | AUC/MIC=%.1f | TCmin=%.1f h | CFUΔ=%.3f\n', ...
            row.PtID, j, dose_mg_per_kg, target.MIC, Cmax, Cmax/target.MIC, AUC/target.MIC, TCmin, CFUBact(j));


        % ---- safety exits so it can NEVER run forever ----
        if j >= maxIters || total_hours_therapy > maxHours
            warning('Safety stop: j=%d, total_hours_therapy=%.1f h (decline never reached 3).', ...
                    j, total_hours_therapy);
            fileID = fopen('./SavedOptimization/Latest/Optimization_failures.txt','a');
            fprintf(fileID, '[Pt %d]\t[Iter %d] Total hours of therapy %3d | CFUΔ=%.3f\n', double(row.PtID), j, total_hours_therapy,CFUBact(j));
            fclose(fileID);
            break
        end

        if total_hours_therapy > 72 && row.cooling == 1, row.cooling = 0; end

        common = intersect(resultOpt.Properties.VariableNames, row.Properties.VariableNames, 'stable');
        resultOpt(j, common) = row(1,common);

        resultOpt.MIC(j)             = target.MIC;
        resultOpt.GFR_L_per_min(j)   = CL_tot_mL_min/1000;
        resultOpt.DoseOpt_mg_kg(j)   = dose_mg_per_kg;
        resultOpt.Dose_mg(j)         = dose_mg;
        resultOpt.Cmax(j)            = Cmax;
        resultOpt.TCmin(j)           = TCmin;
        resultOpt.CmaxMIC(j)         = Cmax/target.MIC;
        resultOpt.AUCMIC(j)          = AUC/target.MIC;
        resultOpt.CFUBact(j)         = CFUBact(j);

        %row.PNA = row.PNA + round(TCmin/24);
        updated_row = updateVirtualPatient(row, TCmin);
        row = updated_row;
        close all
    end

    resultsOptPop = [resultsOptPop; resultOpt(~ismissing(resultOpt.PtID), :)];
end

save_path = './SavedOptimization/Latest';
save(fullfile(save_path, 'Opt_MultiDosePop_AMICREA_CI.mat'), 'resultsOptPop');
writetable(resultsOptPop,fullfile(save_path,"Opt_MultiDosePop_AMICREA_CI.xlsx"));

%% Scheduling/PD results plot
T = resultsOptPop;                          
if ~iscategorical(T.PtID)
    T.PtID = categorical(T.PtID);
end

[G, pid] = findgroups(T.PtID);
T.Cycle  = zeros(height(T),1);
for k = 1:max(G)
    idx = find(G==k);
    T.Cycle(idx) = 1:numel(idx);
end

% ----- Map max cycle number back to each row -----
lastCyclePerGroup = splitapply(@max, T.Cycle, G);   % one value per group
lastCycleForRow   = lastCyclePerGroup(G);           % expand to each row
isLast = T.Cycle == lastCycleForRow;                % TRUE if this row is the last for its patient

% ----- Baseline per patient and relative change -----
baselineCFU = splitapply(@(x) x(1), T.CFUBact, G);  % first value per patient
T.CFUBaseline = baselineCFU(G);
T.RelChange   = T.CFUBact - T.CFUBaseline;

Tfinal = T(isLast,:);
Tfinal.Direction = categorical(Tfinal.RelChange <= 0, [true false], {'Decrease','Increase'});


% (2) Cohort heatmap: patient × cycle of Δlog10 ---------------
T = resultsOptPop;                            % your table with NsLog, NrLog
if ~iscategorical(T.PtID), T.PtID = categorical(T.PtID); end

[G, pid] = findgroups(T.PtID);
T.Cycle  = zeros(height(T),1);
for k = 1:max(G)
    idx = find(G==k);
    T.Cycle(idx) = 1:numel(idx);
end

% Baselines (first row per patient)
baseTot = splitapply(@(x) x(1), T.CFUBact, G);
baseNs  = splitapply(@(x) x(1), T.NsLog,  G);
baseNr  = splitapply(@(x) x(1), T.NrLog,  G);

% Map baselines back to each row
T.CFU0   = baseTot(G);
T.Ns0    = baseNs(G);
T.Nr0    = baseNr(G);

% Relative changes (Δlog10)
T.dTot = T.CFUBact - T.CFU0;   % Total Δlog10
T.dNs  = T.NsLog   - T.Ns0;    % Ns    Δlog10
T.dNr  = T.NrLog   - T.Nr0;    % Nr    Δlog10

% Patient×Cycle matrices (pad ragged with NaN)
[pid,~,pidIdx] = unique(T.PtID,'stable');
maxC = max(T.Cycle);
Mtot = nan(numel(pid), maxC);
Mns  = nan(numel(pid), maxC);
Mnr  = nan(numel(pid), maxC);
for k = 1:numel(pid)
    rows = find(pidIdx==k);
    cyc  = T.Cycle(rows);
    Mtot(k,cyc) = T.dTot(rows);
    Mns(k,cyc)  = T.dNs(rows);
    Mnr(k,cyc)  = T.dNr(rows);
end

% Diverging blue-white-red colormap (centered at 0)
% --- Blue → Light Grey → Red (zero = light grey) ---
Ncol  = 256;
half  = Ncol/2;
blue = [84 39 136]/255;
grey0= [0.92 0.92 0.92];
red = [253 174 97]/255;

cmap1 = [linspace(blue(1),grey0(1),half)', ...
         linspace(blue(2),grey0(2),half)', ...
         linspace(blue(3),grey0(3),half)'];

cmap2 = [linspace(grey0(1),red(1),half)', ...
         linspace(grey0(2),red(2),half)', ...
         linspace(grey0(3),red(3),half)'];

cmap = [cmap1; cmap2];
colormap(cmap)


% ------------------ Plot: 3 panels side-by-side --------------------
fig = figure('Color','w','Position',[120 200 1200 400]);
tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');

grayBG = [1 1 1];
yticksToShow = 1:100:numel(pid);

% ----- Panel 1 : Total -----
ax = nexttile;
h = imagesc(Mtot);                 % no clims or AlphaData here
caxis(clims);
set(h,'AlphaData',~isnan(Mtot));
ax.Color = grayBG;
axis tight; colormap(cmap); colorbar
title('$\Delta \log_{10}$ Total Bacteria ($N_s+N_r$)','Interpreter','latex')
xlabel('Cycle','Interpreter','latex'); ylabel('Patient','Interpreter','latex');
set(ax,'YTick',yticksToShow,'YTickLabel',cellstr(pid(yticksToShow)), ...
       'TickDir','out','Box','on','FontName','Helvetica','FontSize',11);

% ----- Panel 2 : Ns -----
ax = nexttile;
h = imagesc(Mns);
caxis(clims);
set(h,'AlphaData',~isnan(Mns));
ax.Color = grayBG;
axis tight; colormap(cmap); colorbar
title('$\Delta \log_{10}$ Sensitive Bacteria ($N_s$)','Interpreter','latex')
xlabel('Cycle','Interpreter','latex'); ylabel('');
set(ax,'YTick',yticksToShow,'YTickLabel',[], ...
       'TickDir','out','Box','on','FontName','Helvetica','FontSize',11);

% ----- Panel 3 : Nr -----
ax = nexttile;
h = imagesc(Mnr);
caxis(clims);
set(h,'AlphaData',~isnan(Mnr));
ax.Color = grayBG;
axis tight; colormap(cmap);
cb = colorbar;
cb.Label.String = '$\Delta \log_{10}$ CFU relative to baseline';
cb.Label.Interpreter = 'latex';
title('$\Delta \log_{10}$ Resistant Bacteria ($N_r$)','Interpreter','latex')
xlabel('Cycle','Interpreter','latex'); ylabel('');
set(ax,'YTick',yticksToShow,'YTickLabel',[], ...
       'TickDir','out','Box','on','FontName','Helvetica','FontSize',11);


%% ===== Dominance heatmap (absolute logs) =====
Delta0 = (baseNs(1) - baseNr(1));     % > 0 because Ns0 > Nr0

% --- Dominance heatmap on ABSOLUTE logs (zero = equal abundance) ---
% logNr_end - logNs_end = (Mnr - Mns) + (logNr0 - logNs0) = (Mnr - Mns) - Delta0
Mdom = (Mnr - Mns) - Delta0;
yticksToShow = 1:100:numel(pid);      % every 100th patient

% color limits
absVals = abs(Mdom(:)); absVals = absVals(~isnan(absVals));
vmaxD   = isempty(absVals) * 1 + ~isempty(absVals) * max(0.5, prctile(absVals,95));
climsD  = [-vmaxD, vmaxD];

figure('Color','w','Position',[120 200 900 420]);
imagesc(Mdom, climsD);
axis tight

% --- Dark Green → White → Purple colormap ---
n = 256;
darkgreen = [0    143 57]/255;    % darker, richer green
purple     = [255   0   0]/255;
white    = [255 172   28]/255; % vivid purple

half = n/2;
cmap1 = [linspace(darkgreen(1),white(1),half)', ...
         linspace(darkgreen(2),white(2),half)', ...
         linspace(darkgreen(3),white(3),half)'];
cmap2 = [linspace(white(1),purple(1),half)', ...
         linspace(white(2),purple(2),half)', ...
         linspace(white(3),purple(3),half)'];
cmap = [cmap1; cmap2];

colormap(cmap)
caxis(climsD)
colorbar

title('Dominance Heatmap','Interpreter','latex')
xlabel('Cycle','Interpreter','latex')
ylabel('Patient','Interpreter','latex')
set(gca,'YTick',yticksToShow, ...
        'YTickLabel',cellstr(pid_r(yticksToShow)), ...
        'TickDir','out','Box','on','FontName','Helvetica','FontSize',11);

cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10} N_r - \log_{10} N_s$';
cb.Label.FontSize=10;
% Put explanatory text on the top and bottom ticks
cb.Ticks = [climsD(1) 0 climsD(2)];
cb.TickLabelInterpreter = 'latex';
cb.TickLabels = { '$N_s > N_r$', '$N_s=N_r$', '$N_r > N_s$'};
% Overlay light grid
hold on
for c = 1.5:nC-0.5, plot([c c],[0.5 nP+0.5],':','Color',[0 0 0 0.08]); end
for r = 1.5:nP-0.5, plot([0.5 nC+0.5],[r r],':','Color',[0 0 0 0.08]); end

% Hatch/gray NaNs
[nr, nc] = size(Mdom_r);
[nanr, nanc] = find(isnan(Mdom_r));
for k = 1:numel(nanr)
   rectangle('Position',[nanc(k)-0.5, nanr(k)-0.5, 1, 1], ...
             'FaceColor',[1 1 1], 'EdgeColor','none', 'Curvature',0);
end
uistack(findobj(gca,'Type','Image'),'bottom'); % keep image behind overlays
hold off


%% helpers
function val = getparamvalue(mObj, name, defaultIfMissing)
    p = sbioselect(mObj,'Type','parameter','Name',name);
    if isempty(p), val = defaultIfMissing; else, val = p.Value; end
end

function setparamLh(mObj,name,valLh)
    p = sbioselect(mObj,'Type','parameter','Name',name);
    if isempty(p)
        addparameter(mObj,name,valLh,'ValueUnits','liter/hour');
    else
        p.Value = valLh;
        try, p.ValueUnits = 'liter/hour'; end
    end
end

function setparam(mObj,name,val)
    p = sbioselect(mObj,'Type','parameter','Name',name);
    if isempty(p)
        addparameter(mObj,name,val,'ValueUnits','dimensionless');
    else
        p.Value = val;
    end
end


%% helper functions
function [Xopt, fval] = optimizeRegimen(mObj, cs, target, cfg, BW_kg, start)

    % Define bounds for continuous variables
    dose_bounds = cfg.bounds.DoseAmount_mg_per_kg;

    bestSSE = Inf;
    bestX = NaN(1);
    opts = optimoptions('fmincon', 'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 300, 'UseParallel',true,...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-4, ...
    'OptimalityTolerance', 1e-3 ...
    );


    doseInit     = dose_bounds(1)     + rand * diff(dose_bounds);
    X0 = doseInit;
    lb = dose_bounds(1);
    ub = dose_bounds(2);

    objWrap = @(x) objFun(x, start);
    conWrap = @(x) nonlincon(x, start);

    [xOptCont, fval_tmp] = fmincon(objWrap, X0, [], [], [], [], lb, ub, conWrap, opts);
    if fval_tmp < bestSSE
        bestSSE = fval_tmp;
        bestX = [xOptCont(1)];
    end


    Xopt = bestX;
    fval = bestSSE;

    function sse = objFun(X, start)
        [AUC, Cmax, TCmin, ~, ~] = simulateMetrics(mObj, cs, X, cfg, BW_kg, start);
        fu=0.9;
        ratio = fu * Cmax / target.MIC;
    
        % Penalità soft
        e_tcmin = max(0, TCmin - target.TCmin_max);      % TCmin > 25
        %e_cmax  = max(0, ratio  - target.CmaxMIC_lb);  % Cmax > 35
    
        % (opzionale) micro-penalità per ratio fuori banda
        e_ratio_lo = max(0, target.CmaxMIC_lb - ratio);
        % e_ratio_hi = max(0, ratio - target.CmaxMIC_ub);
    
        w = struct('tc',1,'cmax',1);  % pesi principali
        sse = w.tc*e_tcmin^2 + w.cmax*e_ratio_lo^2;
        % sse = sse + 0.5*(e_ratio_lo^2 + e_ratio_hi^2);  % opzionale
    end



    % --- Nested nonlinear constraints ---
    function [c, ceq] = nonlincon(X, start)
        [AUC, Cmax, TCmin, ~, ~] = simulateMetrics(mObj, cs, X, cfg, BW_kg, start);
        ratio = Cmax / target.MIC;
    
        % disuguaglianze c<=0:
        c = [
            %TCmin - target.TCmin_max            % TCmin ≤ 25
            Cmax  - target.Cmax_cap_hard        % Cmax ≤ 60
            %target.CmaxMIC_lb - ratio           % ratio ≥ 8
            % ratio - target.CmaxMIC_ub           % ratio ≤ 10
        ];
        ceq = [];
    end


end



function [AUC, Cmax, TCmin, next_initial_conditions,sd] = simulateMetrics(mObj, cs, X, cfg, BW_kg, initial_conditions)
    
    delete(getdose(mObj, 'Amikacin'));

    % Add new dose
    doseObj = sbiodose('Amikacin', 'schedule');
    doseObj.TargetName  = cfg.doseName;
    dose_mg = X(1) * BW_kg;
    doseObj.Amount = dose_mg; 
    doseObj.AmountUnits = 'milligram';
    doseObj.Rate = dose_mg / cfg.infDur;
    doseObj.RateUnits   = cfg.rateUnits;
    doseObj.Time   = 0;
    doseObj.TimeUnits = 'hour';
    doseObj.Active = true;
    
    % Simulation
    cs.StopTime = cfg.postWindow;
    cs.TimeUnits = 'hour';
    cs.SolverOptions.OutputTimes = 0:0.1:cs.StopTime;  % ogni 6 minuti

    cs.RuntimeOptions.StatesToLog = {'Serum_blood.Cb','Serum_blood.Ns','Serum_blood.Nr','u'};

    if ~isempty(initial_conditions)
        mObj.Species(1,1).InitialAmount = initial_conditions(1); 
        mObj.Species(3,1).InitialAmount = initial_conditions(2); 
        mObj.Species(4,1).InitialAmount = initial_conditions(3);
    else
        mObj.Species(3,1).InitialAmount = 8.8087e+05;
        mObj.Species(4,1).InitialAmount = 2.2194e+05;
    end
    
    cfu_initial_count = log10(8.8087e+05 + 2.2194e+05);

    % Run simulation on time-variant model
    cs.SolverOptions.AbsoluteToleranceScaling = false;
    success = false;

    while ~success
        try
            sd = sbiosimulate(mObj, cs, doseObj);
            success = true; 
        catch ME
            if contains(ME.message, 'step size')
                fprintf('Riduco Post Window e riprovo...\n');
                cfgObj = mObj.addconfigset(string(rand()));
                cfgObj.StopTime = 20; 
                cfgObj.SolverType = 'ode45';
                sd = sbiosimulate(mObj, cfgObj, doseObj);
                cfgObj.delete;
                success = true;
            end
        end
    end

   % Compute PK metrics
    T = sd.Time; % absolute time
    C = sd.Data(:, strcmp(sd.DataNames, cfg.concName));
    [Cmax,Tmax] = max(C);
    tp = find(C(Tmax : end)<=3);

    retry = 0;
    while isempty(tp)
        retry = retry+1 ;
        figure(1);plot(sd.Time, sd.Data(:,1));hold on; yline(3); 
        cfgObj = mObj.addconfigset(string(rand()));
        cfgObj.StopTime = cs.StopTime+retry*10; 
        cfgObj.SolverType = 'ode45';
        sd = sbiosimulate(mObj, cfgObj, doseObj);
        cfgObj.delete;
        % Compute PK metrics
        T = sd.Time; % absolute time
        C = sd.Data(:, strcmp(sd.DataNames, cfg.concName));
        [Cmax,Tmax] = max(C);
        tp = find(C(Tmax : end)<=3);
    end

    tpCmin = Tmax+tp(1);
    TCmin = T(tpCmin);
    AUC   = trapz(T(1:tpCmin), C(1:tpCmin));   % AUC over interdose interval ONLY


    % Final conditions
    next_initial_conditions(1) = sd.Data(tpCmin,1);
    next_initial_conditions(2) = sd.Data(tpCmin,2);
    next_initial_conditions(3) = sd.Data(tpCmin,3);
    % Observe CFU count with respect to starting condition
    next_initial_conditions(4) = (log10(sd.Data(tpCmin,2)+sd.Data(tpCmin,3)) - cfu_initial_count)/cfu_initial_count;
    % next_initial_conditions(5) = sbioselect(mObj, 'Name','u','Type','parameter').Value;
    next_initial_conditions(5) = sd.Data(tpCmin,4);

end

function mic = MIC_from_u(u,baselineMIC)
    u      = max(u,0);
    u50        = 300;     % tune here
    h          = 1;
    frac   = (u.^h) ./ (u50.^h + u.^h);     % in [0,1)
    mic    = baselineMIC + (8 - baselineMIC) * frac;
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

function [bestTau, bestAlpha] = calibrate_alpha_tau(GA, PNA, y, mu_fun, tau_min, tau_max, tau_step)
% Robust 1D search over tau with closed-form alpha per tau in log-space.
    y = double(y); mask = isfinite(y) & (y>0) & isfinite(GA) & isfinite(PNA);
    GA = GA(mask); PNA = PNA(mask); y = y(mask);
    logy = log(y);

    taus = tau_min:tau_step:tau_max;
    bestSSE = Inf; bestTau = 0; bestAlpha = 0;

    for t = taus
        mu = mu_fun(GA, PNA + t);
        mu = max(mu, eps);
        r  = logy - log(mu);
        alpha_hat = mean(r);                 % closed-form minimizer
        SSE = sum( (r - alpha_hat).^2 );     % residual SSE
        if SSE < bestSSE
            bestSSE = SSE; bestTau = t; bestAlpha = alpha_hat;
        end
    end
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

        Cr = min(max(Cr_tot, 0), 2.5);
        yhat(i) = Cr;
    end
end

function v = mpe(yhat,y),  v = mean(yhat - y, 'omitnan'); end
function v = mape(yhat,y), v = mean(abs((yhat - y)./y), 'omitnan')*100; end
function v = rmse_func(yhat,y), v = sqrt(mean((yhat - y).^2, 'omitnan')); end

function T = stratify_bins(x, yhat, y, edges, label)
    [g,~]=discretize(x, edges); K = numel(edges)-1;
    Bin = strings(K,1); N = zeros(K,1); MPE=nan(K,1); MAPE=nan(K,1); RMSE=nan(K,1);
    for k=1:K
        idx = (g==k) & isfinite(y) & isfinite(yhat) & y>0;
        N(k)=sum(idx);
        if N(k)>0
            MPE(k)  = mean(yhat(idx)-y(idx));
            MAPE(k) = mean(abs((yhat(idx)-y(idx))./y(idx)))*100;
            RMSE(k) = sqrt(mean((yhat(idx)-y(idx)).^2));
        end
        Bin(k)=sprintf('%s %g–%g', label, edges(k), edges(k+1));
    end
    T = table(Bin, N, MPE, MAPE, RMSE);
end

function [bestTau, bestAlpha] = calibrate_alpha_tau_pen(GA, PNA, y, mu_fun, tau_min, tau_max, tau_step, lambda)
% Like calibrate_alpha_tau, but adds ridge penalty lambda*tau^2 in log-space SSE.
    y = double(y);
    m = isfinite(y) & y>0 & isfinite(GA) & isfinite(PNA);
    GA = GA(m); PNA = PNA(m); y = y(m); logy = log(y);

    taus = tau_min:tau_step:tau_max;
    bestSSE = Inf; bestTau=0; bestAlpha=0;
    for t = taus
        mu = mu_fun(GA, PNA + t); mu = max(mu, eps);
        r  = logy - log(mu);
        alpha_hat = mean(r);
        SSE = sum((r - alpha_hat).^2) + lambda*(t^2);  % penalize tau^2
        if SSE < bestSSE
            bestSSE = SSE; bestTau=t; bestAlpha=alpha_hat;
        end
    end
end

function cv = kfold_calibrate_eval(K, GA, PNA, y, mu_fun, tau_min, tau_max, tau_step, lambda)
% K-fold CV for calibration: fit (alpha,tau) on train, evaluate on test.
    n = numel(y); idx = find(isfinite(y) & y>0 & isfinite(GA) & isfinite(PNA));
    idx = idx(randperm(numel(idx)));
    folds = ceil((1:numel(idx)) * K/numel(idx));
    cv.paper.MPE=zeros(K,1); cv.paper.MAPE=zeros(K,1); cv.paper.RMSE=zeros(K,1);
    cv.cal.MPE  =zeros(K,1); cv.cal.MAPE  =zeros(K,1); cv.cal.RMSE  =zeros(K,1);

    for k=1:K
        test  = idx(folds==k);
        train = idx(folds~=k);
        % fit on train (penalized)
        [tau_k, alpha_k] = calibrate_alpha_tau_pen(GA(train), PNA(train), y(train), mu_fun, tau_min, tau_max, tau_step, lambda);

        % evaluate on test
        yhat_p = mu_fun(GA(test), PNA(test));
        yhat_c = exp(alpha_k) .* mu_fun(GA(test), PNA(test)+tau_k);
        yo     = y(test);

        cv.paper.MPE(k)  = mpe(yhat_p, yo);
        cv.paper.MAPE(k) = mape(yhat_p, yo);
        cv.paper.RMSE(k) = rmse_func(yhat_p, yo);

        cv.cal.MPE(k)  = mpe(yhat_c, yo);
        cv.cal.MAPE(k) = mape(yhat_c, yo);
        cv.cal.RMSE(k) = rmse_func(yhat_c, yo);
    end
end
function row = band_row(label, idx, GA, PNA, y_obs, paper_mu, B)
    if nnz(idx)==0
        row = table(string(label), 0, NaN, NaN, NaN, NaN, NaN, ...
            'VariableNames', {'Band','N','Med_obs','Med_paper','Ratio','Ratio_lo95','Ratio_hi95'});
        return;
    end
    yo  = y_obs(idx);
    mu  = paper_mu(GA(idx), PNA(idx));
    mu  = max(mu, eps);

    med_obs   = median(yo, 'omitnan');
    med_paper = median(mu, 'omitnan');
    ratio     = med_obs / med_paper;

    % bootstrap CI for ratio of medians
    n = numel(yo);
    rBoot = nan(B,1);
    for b=1:B
        ii = randi(n, n, 1);
        obs_b = median(yo(ii), 'omitnan');
        pap_b = median(mu(ii), 'omitnan');
        rBoot(b) = obs_b / max(pap_b, eps);
    end
    ci = quantile(rBoot, [0.025 0.975]);

    row = table(string(label), nnz(idx), med_obs, med_paper, ratio, ci(1), ci(2), ...
        'VariableNames', {'Band','N','Med_obs','Med_paper','Ratio','Ratio_lo95','Ratio_hi95'});
end
function v = mpe_mask(yhat,y,mask)
    v = mean(yhat(mask) - y(mask), 'omitnan');
end

function v = mape_mask(yhat,y,mask)
    ym = y(mask);
    v  = mean( abs((yhat(mask) - ym)./ym), 'omitnan') * 100;
end

function v = rmse_mask(yhat,y,mask)
    v = sqrt(mean( (yhat(mask) - y(mask)).^2, 'omitnan'));
end

function T = stratify_bins_mask(x, yhat, y, edges, label, mask)
    % Apply common mask first
    x = x(mask); yhat = yhat(mask); y = y(mask);
    [g,~] = discretize(x, edges);
    K = numel(edges)-1;
    Bin = strings(K,1); N = zeros(K,1); MPE=nan(K,1); MAPE=nan(K,1); RMSE=nan(K,1);
    for k=1:K
        idx = (g==k) & isfinite(y) & isfinite(yhat) & (y>0);
        N(k) = sum(idx);
        if N(k)>0
            dy = yhat(idx) - y(idx);
            MPE(k)  = mean(dy);
            MAPE(k) = mean(abs(dy ./ y(idx))) * 100;
            RMSE(k) = sqrt(mean(dy.^2));
        end
        Bin(k) = sprintf('%s %g–%g', label, edges(k), edges(k+1));
    end
    T = table(Bin, N, MPE, MAPE, RMSE);
end

function M = metrics_bundle_mask(yhat, y, mask)
% Full bundle: multiplicative + absolute + bias + legacy MAPE.
    yh = double(yhat(mask)); yo = double(y(mask));
    % Guard
    m = isfinite(yh) & isfinite(yo) & yh>0 & yo>0;
    yh = yh(m); yo = yo(m);
    N  = numel(yh);

    % log-space
    loge  = log(yh) - log(yo);
    RMSLE = sqrt(mean(loge.^2));
    MALE  = mean(abs(loge));
    GMAPE = (exp(MALE) - 1) * 100;           % %
    GMR   = exp(mean(loge));                 % ratio yhat/y

    % Bootstrap CI for GMR
    B = 1000;
    if exist('bootstrp','file')
        a = bootstrp(B, @mean, loge);
    else
        a = zeros(B,1);
        n = numel(loge);
        for b=1:B
            idx = randi(n, n, 1);
            a(b) = mean(loge(idx));
        end
    end
    ci = quantile(a, [0.025 0.975]);
    GMR_lo95 = exp(ci(1)); GMR_hi95 = exp(ci(2));

    % absolute (mg/dL) + robust
    abs_err = abs(yh - yo);
    MAE  = mean(abs_err);
    RMSE = sqrt(mean((yh - yo).^2));
    MedAE= median(abs_err);

    % sMAPE (%)
    sMAPE = mean( 2*abs(yh-yo) ./ max(abs(yh)+abs(yo), eps) ) * 100;

    % legacy MPE/MAPE for continuity
    MPE  = mean(yh - yo);
    MAPE = mean(abs((yh - yo)./yo)) * 100;

    M = struct('N',N,'RMSLE',RMSLE,'MALE',MALE,'GMAPE_pct',GMAPE, ...
               'GMR',GMR,'GMR_lo95',GMR_lo95,'GMR_hi95',GMR_hi95, ...
               'MAE',MAE,'RMSE',RMSE,'MedAE',MedAE,'sMAPE_pct',sMAPE, ...
               'MPE',MPE,'MAPE',MAPE);
end

function C = calibration_log_mask(yhat, y, mask)
% Regress log(y) on log(yhat): intercept≈0, slope≈1 if calibrated.
    yh = double(yhat(mask)); yo = double(y(mask));
    m  = isfinite(yh) & isfinite(yo) & yh>0 & yo>0;
    X  = log(yh(m)); Y = log(yo(m));
    if isempty(X)
        C = struct('intercept',NaN,'slope',NaN,'R2',NaN);
        return;
    end
    mdl = fitlm(X, Y, 'Intercept', true, 'RobustOpts', 'on');
    C = struct('intercept', mdl.Coefficients.Estimate(1), ...
               'slope',     mdl.Coefficients.Estimate(2), ...
               'R2',        mdl.Rsquared.Ordinary);
end


function [tout, yout] = simulate_pbpk_one(p, BW_kg, dose_mgkg, dose_times_h, t_end_h, y0, theta)
    % dose (mg) into depot; depot empties with p.kf ≈ 1-h infusion
    p.EC50s0 = theta(1); p.EC50r0 = theta(2); p.c = theta(3);
    p.rS     = theta(4); p.Emax      = theta(5); p.rR = theta(6); p.ku=theta(7); p.sigma=theta(8); p.d=theta(9);

    % Fit ['EC50s0','EC50r0','c','rS','Emax','rR','ku','sigma','d'] for a *fixed* physiology p



    % keep all scheduled doses within horizon (including 0 h)
    dose_times = unique(dose_times_h(dose_times_h <= t_end_h));
    % positive (nonzero) doses for segmentation
    pos_doses  = dose_times(dose_times > 0);

    dose_mg    = dose_mgkg * BW_kg;

    opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'NonNegative',1:8);
    tout = []; yout = [];
    y = y0;

    % give dose at t=0 if scheduled
    if any(abs(dose_times - 0) < 1e-12)
        y(2) = y(2) + dose_mg;
    end

    % integrate from segment to segment; add dose right after each pos_dose time
    all_times  = [pos_doses, t_end_h];
    seg_starts = [0,        pos_doses];

    for i = 1:numel(all_times)
        tspan = [seg_starts(i), all_times(i)];
        sol = ode45(@(t,yy) rhs_pbpk_pd(t,yy,p), tspan, y, opts);
        tt  = linspace(tspan(1), tspan(2), max(2, ceil((tspan(2)-tspan(1))*6)));
        yy  = deval(sol, tt).';
        tout = [tout; tt(:)];
        yout = [yout; yy];

        % if this segment ended at a scheduled (positive) dose time, give the dose
        if i <= numel(pos_doses)
            y = yy(end,:).';
            y(2) = y(2) + dose_mg;      % add next dose into depot
        else
            y = yy(end,:).';
        end
    end

    [tout, idx] = unique(tout,'stable'); 
    yout = yout(idx,:);
end


function dy = rhs_pbpk_pd(~, y, p)
Cb=y(1); Amino=y(2); Ch=y(3); Ck=y(4); Cr=y(5); u=y(6); Ns=y(7); Nr=y(8);
Sb=p.Serum_blood; He=p.Hearth; Ki=p.Kidney; Re=p.Remainder;
Qh=p.BloodFlowRateHearth; Qk=p.BloodFlowRateKidney; Qr=p.BloodFlowRateRemainder;
Ph=p.PartitionCoefficientHearth; Pk=p.PartitionCoefficientKidney; Pr=p.PartitionCoefficientRemainder;
CL=p.CL; kf=p.kf;
EC50s0=p.EC50s0; EC50r0=p.EC50r0; c=p.c; rS=p.rS; Emax=p.Emax; rR=p.rR;
sigma=p.sigma; ku=p.ku; aRS=p.aRS; aSR=p.aSR; Kcap=p.K; d=p.d;
expInput = ku*u; Cf = p.fu * Cb;

dCb = (1/Sb) * ( ...
      -(Qk*Cb)*Sb + (Qk*(Ck/Pk))*Ki ...
    -  (Qh*Cb)   +  (Qh*(Ch/Ph)) ...
    +  (Qr*(Cr/Pr))*Re - (Qr*Cb)*Sb ...
    +   (kf*Amino) - (CL*Cb) );
dAmino = -(kf*Amino);
dCh = (1/He) * ( Qh*Cb - Qh*(Ch/Ph) );
dCk = (1/Ki) * ( (Qk*Cb)*Sb - (Qk*(Ck/Pk))*Ki );
dCr = (1/Re) * ( -(Qr*(Cr/Pr))*Re + (Qr*Cb)*Sb );
du  = sigma * ( - rR*c*exp(-c*u)*(1 - (Nr + aRS*Ns)/Kcap) ...
      + Emax * Cf * EC50r0 * ku * exp(expInput) / ( (EC50r0*exp(expInput) + Cf)^2 ) );
dNs = ( rS*(1 - (Ns + aSR*Nr)/Kcap) - Emax*Cf/(EC50s0 + Cf) - d ) * Ns;
dNr = ( rR*exp(-c*u)*(1 - (Nr + aRS*Ns)/Kcap) - Emax*Cf/(EC50r0*exp(expInput) + Cf) - d ) * Nr;

dy = [dCb; dAmino; dCh; dCk; dCr; du; dNs; dNr];
end


function GAflag = ga_bucket_from_row(row)
% 25 for GA<32, 40 otherwise (matches your description)
    GAflag = 25;
    if isfield(row,'GA')
        if row.GA >= 32, GAflag = 40; end
    elseif isprop(row,'GA')
        if row.GA >= 32, GAflag = 40; end
    else
        error('Row has no GA field.');
    end
end

function M = genta_csv_metrics(baseDir, sched, GAflag)
% sched: 'g1'|'g2'|'g3'
% GAflag: 25 or 40
% CSV format: 2 cols [time_h, conc_mgL]
    dose_map = struct('g1',4,   'g2',4.5, 'g3',5);
    tau_map  = struct('g1',24,  'g2',36,  'g3',48);

    fn = fullfile(baseDir, sprintf('%s_%d.csv', sched, GAflag));
    if ~isfile(fn), error('CSV not found: %s', fn); end
    T = readtable(fn);
    t = double(T{:,1});
    C = double(T{:,2});
    [Cmax, iMax] = max(C);

    % TCmin (first ≤3 mg/L after Cmax). Loosen if never crosses.
    ii = find(C(iMax:end) <= 3, 1, 'first');
    if isempty(ii), TCmin = t(end); else, TCmin = t(iMax+ii-1); end

    M = struct();
    M.file = fn;
    M.sched = sched;
    M.GAflag = GAflag;
    M.t = t(:);
    M.C = C(:);
    M.Cmax = Cmax;
    M.TCmin = TCmin;
    M.dose_mgkg_gent = dose_map.(sched);
    M.tau_h = tau_map.(sched);

    % Reference PD driver for gent; change fu/MIC if you have better values
    MIC_gent = 2; fu_gent = 1.0;
    M.fCmax_MIC_ref = (fu_gent*Cmax)/MIC_gent;
end

function [schedules, tag] = chooseTargetsByGA(GA)
% GA<=32 -> GA=25 targets; GA>32 -> GA=40 targets
if GA <= 32
    schedules = {'g1_25','g2_25','g3_25'}; tag = 'GA25';
else
    schedules = {'g1_40','g2_40','g3_40'}; tag = 'GA40';
end
end

function [Cb0, CbT] = get_Cb_start_end(simData)
    % simData: SimBiology SimData object from sbiosimulate
    % Returns:
    %   Cb0 = Cb concentration at simulation start
    %   CbT = Cb concentration at simulation end

    % time vector
    t  = simData.Time;
    % find which column corresponds to Cb
    id = strcmp(simData.DataNames, 'Cb');
    % pull out that trajectory
    Cb = simData.Data(:, id);

    % first and last values
    Cb0 = Cb(1);
    CbT = Cb(end);
end

function resetInitialsToSimData(modelObj, simData)
    % Carry states forward: set all species/compartment init values to last point
    sps = sbioselect(modelObj,'Type','species');
    for ii = 1:numel(sps)
        nm = sps(ii).Name;
        j  = find(strcmp(simData.DataNames, nm), 1, 'first');
        if ~isempty(j)
            sps(ii).InitialAmount = simData.Data(end, j);
        end
    end
end

function d = makeInfusionDose(~, name, amt_mg, infDur_h, targetName, amtUnits, rateUnits)
    % Return a standalone ScheduleDose (not attached to the model).
    d = sbiodose(name, 'schedule');
    d.TargetName  = targetName;        % e.g., 'Serum_blood.Aminoglycoside'
    d.Amount      = amt_mg;
    d.AmountUnits = amtUnits;          % 'milligram'
    d.Rate        = amt_mg / infDur_h; % mg/hour over the infusion
    d.RateUnits   = rateUnits;         % 'milligram/hour'
    d.Time        = 0;                 % time from sim start
    d.TimeUnits   = 'hour';
    d.Active      = true;
end


function out = firstFiniteOr(t, var, fallback)
% FIRSTFINITEOR Safely return the first finite numeric value found in table column `var`.
% - Works if the column is numeric, cellstr, string, or categorical.
% - Returns `fallback` if the column is missing or has no finite values.

    out = fallback;

    if ~istable(t) || ~ismember(var, t.Properties.VariableNames)
        return;
    end

    x = t.(var);

    % Normalize to numeric
    if isnumeric(x)
        % ok
    elseif iscell(x) || isstring(x) || ischar(x) || iscategorical(x)
        x = str2double(string(x));
    else
        % unsupported type → give up
        return;
    end

    % Keep finite values only
    x = x(isfinite(x));

    if ~isempty(x)
        out = x(1);
    end
end



function Vdpk = pickVdPerKg(GA, PNA)
    if GA <= 28
        if      PNA < 7,   Vdpk = 0.414;
        elseif  PNA <= 28, Vdpk = 0.472;
        else,              Vdpk = 0.353;
        end
    elseif GA <= 33
        if      PNA < 7,   Vdpk = 0.462;
        else,              Vdpk = 0.454;
        end
    else
        Vdpk = 0.433;
    end
end


function m = tryCopyModel(modelObj)
    if exist('sbiocopy','file')
        m = sbiocopy(modelObj);
    else
        m = copyobj(modelObj);   % fallback for older MATLAB
    end
end

function setCapL(m, compName, valL)
    c = sbioselect(m,'Type','compartment','Name',compName);
    if ~isempty(c) && isfinite(valL) && valL>0
        c.Capacity = valL;
        try, c.CapacityUnits = 'liter'; end %#ok<TRYNC>
    end
end

function v = getParamOr(m, name, default)
    p = sbioselect(m,'Type','parameter','Name',name);
    if isempty(p), v = default; else, v = p.Value; end
end

function setParamIfExists(m, name, value)
    p = sbioselect(m,'Type','parameter','Name',name);
    if ~isempty(p) && isfinite(value)
        p.Value = value;
    end
end


function doseTarget = resolveDoseTarget(m, defaultPath)
% Try default "Comp.Species"; otherwise find plausible alternatives.
    doseTarget = "";

    % 1) Try default path "Comp.Spec"
    parts = split(string(defaultPath), '.');
    if numel(parts) == 2
        comp = parts(1); spec = parts(2);
        compObj = sbioselect(m, 'Type','compartment', 'Name', char(comp));
        if ~isempty(compObj)
            s = sbioselect(compObj, 'Type','species', 'Name', char(spec));
            if ~isempty(s)
                doseTarget = string(comp) + "." + string(spec);
                return;
            end
        end
    end

    % 2) Try common alternatives (Comp.Spec)
    cand = [ "Serum_blood.Amikacin", "Serum_blood.Aminoglycoside", ...
             "Central.Amikacin", "Central.Aminoglycoside", ...
             "Plasma.Amikacin", "Plasma.Aminoglycoside" ];
    for k = 1:numel(cand)
        parts = split(cand(k), '.');
        if numel(parts) ~= 2, continue; end
        comp = parts(1); spec = parts(2);
        compObj = sbioselect(m, 'Type','compartment', 'Name', char(comp));
        if isempty(compObj), continue; end
        s = sbioselect(compObj, 'Type','species', 'Name', char(spec));
        if ~isempty(s)
            doseTarget = cand(k);
            return;
        end
    end

    % 3) Heuristic: any species whose full path contains "ami"
    allSpecies = sbioselect(m, 'Type','species');
    for k = 1:numel(allSpecies)
        full = string(allSpecies(k).Parent.Name) + "." + string(allSpecies(k).Name);
        if contains(lower(full), "ami")
            doseTarget = full;
            return;
        end
    end
end



function yq = safeInterp(tSim, ySim, tObs)
    % Robust interpolation that tolerates duplicates/NaNs/unordered times.
    tSim = double(tSim(:));
    ySim = double(ySim(:));
    tObs = double(tObs(:));

    good = isfinite(tSim) & isfinite(ySim);
    tSim = tSim(good);
    ySim = ySim(good);

    if isempty(tSim)
        yq = nan(size(tObs));
        return;
    end

    [tSim, ord] = sort(tSim);
    ySim = ySim(ord);

    % keep the last value at each duplicated time (or use @mean if you prefer)
    [tU, idxLast] = unique(tSim, 'last');
    yU = ySim(idxLast);

    if numel(tU) == 1
        yq = repmat(yU, size(tObs));
        return;
    end

    yq = interp1(tU, yU, tObs, 'linear', 'extrap');
end

