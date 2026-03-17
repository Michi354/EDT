function out = computePooledFoldMetrics(dataset, modelObj, useTV)
% Pooled fold-error metrics (AAFE headline) with optional time-varying physiology.
% Mirrors evaluateAmikacinAccuracy logic so pooled results are consistent.

    doseTargetDefault = "Serum_blood.Aminoglycoside";
    respNameDefault   = "Serum_blood.Cb";

    % Normalize PtID
    PtID = dataset.PtID;
    if iscell(PtID),        PtID = string(PtID); end
    if iscategorical(PtID), PtID = string(PtID); end
    IDs = unique(PtID,'stable');

    obsAll  = [];
    predAll = [];

    for i = 1:numel(IDs)
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
            continue;
        end

        % ---- GA/PNA for Vd rules ----
        GA  = firstFiniteOr(sub,'GA',NaN);
        PNA = firstFiniteOr(sub,'PNA',NaN);

        % ---- Target Vd (L) ----
        Vd_per_kg = pickVdPerKg(GA,PNA);
        Vd_target = Vd_per_kg * BW_kg;

        % ---- Capacities (L) ----
        BV_per_kg = (GA<37)*0.106 + (GA>=37)*0.080;   % central blood/plasma
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
        if isempty(tObs), continue; end

        % ===== model clone =====
        m  = tryCopyModel(modelObj);
        cs = getconfigset(m);
        cs.StopTime = max(tObs)+24;

        % Allow events to write to these parameters (if useTV)
        if useTV, unlockForEvents(m); end

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

        % ---- Push baseline flows/CL (no EliminationRate to avoid double-clearance) ----
        setParamIfExists(m,'BloodFlowRateHearth',     Qh);
        setParamIfExists(m,'BloodFlowRateHeart',      Qh);
        setParamIfExists(m,'BloodFlowRateKidney',     Qk);
        setParamIfExists(m,'BloodFlowRateRemainder',  Qr);
        setParamIfExists(m,'CL',                      CL_L_h);
        setParamIfExists(m,'TBW',                     BW_kg);
        setParamIfExists(m,'GFR',                     eGFRkg);

        % ---- Dose target ----
        doseTarget = resolveDoseTarget(m, doseTargetDefault);
        if doseTarget == "", continue; end

        % ---- TV physiology events ----
        tch = [];
        if useTV
            [~, tch] = addTimeVaryingPhysiology(m, sub, BW_kg, Vc);
        end

        % ---- Make solver hit obs & event times ----
        try
            cs.SolverOptions.OutputTimes = unique([0; double(tObs(:)); double(tch(:)); cs.StopTime]);
        catch
            % old MATLAB may not support OutputTimes – safe to ignore
        end

        % ---- Simulate ----
        try
            if ~isempty(dTim)
                d = sbiodose('Amikacin','schedule');
                d.TargetName  = char(doseTarget);
                d.Time        = dTim(:);
                d.Amount      = dAmt(:);
                d.AmountUnits = 'milligram';
                d.Rate        = d.Amount ./ 1;
                d.RateUnits   = 'milligram/hour';
                d.TimeUnits   = 'hour';
                d.Active      = true;
                sd = sbiosimulate(m, d);
            else
                sd = sbiosimulate(m);
            end
        catch
            continue;    % skip problematic subject
        end

        % ---- Response extraction ----
        dataNames = string(sd.DataNames);
        idx = find(dataNames == respNameDefault, 1);
        if isempty(idx), idx = find(dataNames == "Cb", 1); end
        if isempty(idx), continue; end

        tSim = double(sd.Time(:));
        ySim = double(sd.Data(:, idx));
        yHat = safeInterp(tSim, ySim, double(tObs(:)));

        good = isfinite(yHat) & isfinite(yObs) & yHat>0 & yObs>0;
        if any(good)
            obsAll  = [obsAll;  double(yObs(good))];   %#ok<AGROW>
            predAll = [predAll; double(yHat(good))];   %#ok<AGROW>
        end
    end

    if isempty(obsAll)
        out = struct('AFE', NaN, 'AAFE', NaN, 'P2', NaN, 'Npairs', 0);
        return;
    end

    lr = log10(predAll ./ obsAll);
    out.AFE    = 10.^mean(lr);
    out.AAFE   = 10.^mean(abs(lr));
    out.P2     = 100*mean((predAll./obsAll >= 0.5) & (predAll./obsAll <= 2));
    out.Npairs = numel(lr);
end

% --- helper to allow event assignments (only needed when useTV=true) ---
function unlockForEvents(m)
    names = { ...
        'CL','EliminationRate', ...
        'BloodFlowRateHearth','BloodFlowRateHeart', ...
        'BloodFlowRateKidney','BloodFlowRateRemainder', ...
        'PartitionCoefficientRemainder' ...
    };
    for k = 1:numel(names)
        p = sbioselect(m,'Type','parameter','Name',names{k});
        for j = 1:numel(p)
            try, p(j).Constant = false; end
        end
    end
end


% ================= helpers =================

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

function out = firstFiniteOr(t, var, fallback)
    out = fallback;
    if ismember(var, t.Properties.VariableNames)
        x = t.(var); x = x(isfinite(x));
        if ~isempty(x), out = x(1); end
    end
end

function Vdpk = pickVdPerKg(GA, PNA)
    if GA <= 28                    % like PMA ≤29 when PNA is small
        if      PNA < 7
            Vdpk = 0.414;   % their subgroup 1 mean
        elseif  PNA <= 28
            Vdpk = 0.472;   % subgroup 2 mean
        else               
            Vdpk = 0.353;   % subgroup 3 mean (tiny n)
        end

    elseif GA <= 33                % like PMA 30–34
        if      PNA < 7   
            Vdpk = 0.462;   % subgroup 4 mean
        else               
            Vdpk = 0.454;   % subgroup 5 mean
        end

    else                           % ga ≥ 34  (-> PMA ≥35, all PNA in paper)
        Vdpk = 0.433;                    % subgroup 6 mean
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
