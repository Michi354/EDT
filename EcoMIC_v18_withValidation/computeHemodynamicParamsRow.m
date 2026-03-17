function out = computeHemodynamicParamsRow(row, flag, varargin)
% computeHemodynamicParamsRow
% Single-subject hemodynamics consistent with createVirtualPop():
%   HR (TH vs non-TH, GA-aware), SV_index (GA/PNA/TH + inotrope; PDA effect),
%   SV (mL), CO (L/min), CI (L/min/m^2).
%
% REQUIRED:
%   row  : 1x1 table with at least GA (weeks), BW (g), PNA (days)
%   flag : your pipeline flag (0/1/2/3). See cooling defaults below.
%
% Name-Value options:
%   'Ht'        : height in cm (if omitted, calls deriveHeight_VolumeDistribution)
%   'DerivePDA' : logical (default true); synthesize PDA/severity if absent
%   'Seed'      : integer for reproducibility (default 42, [] leaves RNG)
%
% Returns struct with: HR, SV_mL, CO_L_min, CI_L_min_m2, SV_index, BSA, Ht,
%                      PDA, PDA_sev, usedCooling, usedIno, usedIbu

    % ---- options ----
    p = inputParser;
    addParameter(p, 'Ht', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'DerivePDA', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Seed', 42, @(x) isempty(x) || isnumeric(x));
    parse(p, varargin{:});
    Ht_in     = p.Results.Ht;
    derivePDA = p.Results.DerivePDA;
    seed      = p.Results.Seed;

    if ~isempty(seed), rng(seed); end

    % ---- safe getters (no dynamic access if missing) ----
    hasVar = @(v) ismember(v, row.Properties.VariableNames);
    function val = getOr(rowLocal, v, defaultVal)
        if ismember(v, rowLocal.Properties.VariableNames)
            val = rowLocal.(v);
            if iscell(val) && numel(val)==1, val = val{1}; end
        else
            val = defaultVal;
        end
    end

    % ---- basics ----
    GA    = row.GA;               % weeks
    BWg   = row.BW;               % grams
    BW_kg = BWg/1000;
    PNA_d = row.PNA;              % days
    PNA_h = PNA_d*24;

    % Inotrope/ibuprofen: default 0 if missing
    Ino = getOr(row,'Ino',0);
    Ibu = getOr(row,'Ibu',0);

    % cooling by flag: 0→0, 1→0 (unless present), 2→use column (fallback 0)
    switch flag
        case 0
            cooling = 0;
        case 1
            cooling = 1;  % usually absent, default 0
        case 2
            cooling = getOr(row,'cooling',0);  % expected present
        otherwise
            cooling = getOr(row,'cooling',0);
    end

    % ---- Height / BSA ----
    Ht = Ht_in;
    if isempty(Ht)
        try
            if flag==2 && hasVar('Gender')
                [Ht, ~] = deriveHeight_VolumeDistribution(PNA_d, GA, BWg, row.Gender);
            else
                [Ht, ~] = deriveHeight_VolumeDistribution(PNA_d, GA, BWg);
            end
        catch
            error('Height not provided and deriveHeight_VolumeDistribution failed. Pass ''Ht'' or ensure the function is on the path.');
        end
    end
    BSA = 0.024265 .* (BW_kg .^ 0.5378) .* (Ht .^ 0.3964);  % Haycock (m^2)

    % ---- PDA (use if present; optionally derive) ----
    PDA     = getOr(row,'PDA',0);
    PDA_sev = getOr(row,'PDA_sev',0);

    if derivePDA && (~hasVar('PDA') || ~hasVar('PDA_sev'))
        % baseline PDA probability by GA (as in createVirtualPop)
        if     GA < 28, p0 = 0.9;
        elseif GA < 30, p0 = 0.8;
        elseif GA < 32, p0 = 0.6;
        elseif GA < 34, p0 = 0.5;
        elseif GA < 36, p0 = 0.35;
        elseif GA < 37, p0 = 0.20;
        else,            p0 = 0.01;
        end
        p = p0 * exp(-PNA_d/9);
        if Ibu==1, p = p * 0.6; end

        isPDA = rand < p; sev = 0;
        if isPDA
            r = rand;
            if     r < 0.50, sev = 1;
            elseif r < 0.85, sev = 2;
            else,             sev = 3;
            end
            if PNA_d > 3, sev = max(1, sev-1); end
            if Ibu==1 && sev>1, sev = sev-1; end
        end
        PDA     = double(isPDA);
        PDA_sev = sev;
    end

    % ---- constants (from createVirtualPop) ----
    % HR (TH vs non-TH)
    mu12 = 92; mu72 = 97; sdTH = 8; ino_delta_TH = 15;
    a_HR = 202.667; b_HR = -1.9167; sdNT = 12; ino_delta_NT = 5;

    % SV_index
    sv_base_28w   = 22.0;
    sv_slope_perW = -0.10;
    sv_sd         = 1.8;
    sv_minmax     = [15 28];
    sv_early_PNA48 = +1.0;
    sv_TH_term     = -1.0;
    sv_inotrope    = +1.0;

    sev_gain = [0.00, 0.05, 0.10, 0.18];

    % ---- Heart Rate (bpm) ----
    if (cooling==1) && (GA>=36)
        t = min(max(PNA_h,0),72);
        if t <= 12
            muHR = mu12;
        else
            muHR = mu12 + (mu72-mu12)*(t-12)/(72-12);
        end
        if Ino==1, muHR = muHR + ino_delta_TH; end
        HR = round(muHR + sdTH*randn(1));
        HR = min(max(HR, 80), 120);
    else
        muHR = a_HR + b_HR * GA;
        if PNA_d > 5, muHR = muHR - 3; end
        if Ino==1,    muHR = muHR + ino_delta_NT; end
        HR = round(muHR + sdNT*randn(1));
        HR = min(max(HR, 95), 170);
    end

    % ---- SV_index (mL/m^2) ----
    baseSVI = sv_base_28w + sv_slope_perW*(GA - 28);
    if PNA_d <= 2,               baseSVI = baseSVI + sv_early_PNA48; end
    if (cooling==1) && (GA>=36), baseSVI = baseSVI + sv_TH_term;     end
    if Ino==1,                   baseSVI = baseSVI + sv_inotrope;    end

    SV_index = baseSVI + sv_sd*randn(1);

    % PDA effect (week-1, GA<32) with BW-dependent damping & caps
    if (PDA==1) && (PNA_d <= 6) && (GA < 32)
        decay = exp(-PNA_d/1.4);
        bwkg  = BW_kg;
        if      bwkg < 1.0, damp = 0.55;
        elseif  bwkg < 1.5, damp = 0.55 + 0.30*(bwkg - 1.0)/0.5;
        elseif  bwkg < 2.0, damp = 0.85 + 0.15*(bwkg - 1.5)/0.5;
        else,                damp = 1.00;
        end
        base_extra = sev_gain(PDA_sev+1) * decay * damp;
        max_extra  = (bwkg < 1.25) * 0.12 + (bwkg >= 1.25) * 0.22;
        extra      = min(base_extra, max_extra);
        SV_index   = SV_index * (1 + extra);
    end

    % clamp
    SV_index = min(max(SV_index, sv_minmax(1)), 30);

    % ---- SV / CO / CI ----
    SV_mL       = SV_index * BSA;       % mL
    CO_L_min    = (HR * SV_mL) / 1000;  % L/min
    CI_L_min_m2 = CO_L_min / max(BSA, eps);

    % ---- pack output ----
    out = struct();
    out.HR           = HR;
    out.SV_mL        = SV_mL;
    out.CO_L_min     = CO_L_min;
    out.CI_L_min_m2  = CI_L_min_m2;
    out.SV_index     = SV_index;
    out.BSA          = BSA;
    out.Ht           = Ht;
    out.PDA          = PDA;
    out.PDA_sev      = PDA_sev;
    out.usedCooling  = cooling;
    out.usedIno      = Ino;
    out.usedIbu      = Ibu;
end
