function report = evaluateAmikacinAccuracy(dataset, modelObj, useTV)
% useTV = true  -> addTimeVaryingPhysiology (CL/flows/Kp_rem updated over time)
% useTV = false -> keep physiology fixed at baseline (first finite eGFR/CO/etc.)

    PtID = dataset.PtID;
    if iscell(PtID), PtID = string(PtID); end
    if iscategorical(PtID), PtID = string(PtID); end
    [IDs, ~, ~] = unique(PtID, 'stable');

    N = numel(IDs);
    AFE  = nan(N,1);
    AAFE = nan(N,1);
    RMSE = nan(N,1);
    MAE  = nan(N,1);

    doseTargetDefault = "Serum_blood.Aminoglycoside";
    respNameDefault   = "Serum_blood.Cb";

    for i = 1:N
        id  = IDs(i);
        sub = dataset(PtID==id, :);

        % ---- BW (kg) ----
        BW_kg = firstFiniteOr(sub,'BW_kg_',NaN);
        if ~isfinite(BW_kg) || BW_kg<=0
            if ismember('BW',sub.Properties.VariableNames) && isfinite(sub.BW(1))
                BW_kg = sub.BW(1)/1000;
            else
                BW_kg = 3.63;
            end
        end

        % ---- eGFR (mL/min/kg) and CO (L/min) ----
        eGFRkg   = firstFiniteOr(sub,'eGFR',NaN);
        CO_L_min = firstFiniteOr(sub,'CardiacOutput',NaN);
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
                sd = sbiosimulate(m, d);
            else
                sd = sbiosimulate(m);
            end
        catch ME
            fprintf('Pt %s: simulation failed: %s\n', string(id), ME.message);
            continue;
        end

        % ---- Response extraction ----
        dataNames = string(sd.DataNames);
        idx = find(dataNames == respNameDefault, 1);
        if isempty(idx), idx = find(dataNames == "Cb", 1); end
        if isempty(idx)
            fprintf('Pt %s: response not found\n', string(id));
            continue;
        end

        tSim = double(sd.Time(:));
        ySim = double(sd.Data(:, idx));

        % quick guard: if all zeros/NaN, tell us why
        if ~any(isfinite(ySim))
            fprintf('Pt %s: all NaN predictions\n', string(id));
            continue;
        end

        yHat = safeInterp(tSim, ySim, double(tObs(:)));

        good = isfinite(yHat) & isfinite(yObs) & (yObs > 0) & (yHat > 0);
        if ~any(good)
            fprintf('Pt %s: no valid pairs (min yHat=%.3g, max yHat=%.3g)\n', ...
                    string(id), min(yHat), max(yHat));
            continue;
        end

        err      = yHat(good) - yObs(good);
        logRatio = log10(yHat(good)./yObs(good));
        AFE(i)   = 10.^mean(logRatio);
        AAFE(i)  = 10.^mean(abs(logRatio));
        RMSE(i)  = sqrt(mean(err.^2));
        MAE(i)   = mean(abs(err));
    end

    report = table(IDs, AFE, AAFE, RMSE, MAE, ...
        'VariableNames', {'PtID','AFE','AAFE','RMSE','MAE'});
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
