function updated_row = updateVirtualPatient(row, TCmin)

PNA = row.PNA + (TCmin/24);
row.PNA_h = PNA*24;
row.PNA=round(PNA);

GA = row.GA;
PMA = GA + round(row.PNA./7);
Gender = row.Gender;

days = row.PNA;

if Gender == "M" % neonate is Male

    if GA <= 26
        BW = 1.650186e-07*(days^6)- 3.972552e-05*(days^5) + 3.685725e-03*(days^4) ...
                 - 1.667271e-01*(days^3) + 4.012047e00*(days^2) - 3.705484e01*(days)+ 8.887167e02;
        Ht =  8.905395e-04*(days^2)+ 4.139884e-02*(days)+ 3.375000e01;
        Vd =  0.9*BW;
    
    elseif GA <= 29
        BW = 2.59724e-06*(days^6)- 4.23233e-04*(days^5) + 2.69872e-02*(days^4)...
                 - 8.49142e-01*(days^3)+ 1.38614e01*(days^2)- 9.52883e01*(days) + 1.18552e03;
        Ht =  3.14905e-08*(days^3)- 1.05086e-04*(days^2)+ 1.29879e-01*(days)+ 3.52325e01;
        Vd =  0.85*BW;
    
    elseif GA <= 32
        BW =  1.374450e-05*(days^6)- 1.525030e-03*(days^5)+ 6.816827e-02*(days^4)-...
                  1.566683e00*(days^3)+ 1.979696e01*(days^2)- 1.113418e02*(days)+ 1.565761e03;
        Ht =   3.662092e-08*(days^3)- 1.099299e-04*(days^2)+1.254900e-01*(days)+ 3.901344e01;
        Vd =   0.85*BW;
    
    elseif GA <= 36
        BW =  3.097386e-03*(days^4)- 1.981930e-01*(days^3)+ 5.740936e00*(days^2)-...
                  4.577231e01*(days)+ 1.978709e03;
        Ht =   -3.501005e-09*(days^4)+ 3.405542e-06*(days^3)- 1.218821e-03*(days^2)+ 2.462250e-01*(days)+ 4.342661e01;
        Vd =   0.8*BW;
    
    else 
        BW =  1.076382e-11*(days^5)- 5.095864e-08*(days^4) + 9.191582e-05*(days^3)-...
                  7.825887e-02*(days^2)+ 3.746039e01*(days) + 3.579470e03;
        Ht =   2.873236e-08*(days^3)- 8.491953e-05*(days^2) + 9.965299e-02*(days)+ 5.173054e01;
        Vd =   0.75*BW;
    end

    TBW = 0.79 * BW;  % Males: 79% of BW

else
    
    if GA <= 26
        BW = 1.002907e-07*(days^6)- 2.415770e-05*(days^5)+ 2.234812e-03*(days^4)- 1.010066e-01*(days^3)+...
                 2.532719e00*(days^2) - 2.215122e01*(days) + 8.411017e02;
        Ht =  7.419743e-04*(days^2)+ 5.619694e-02*(days) + 3.343850e01;
        Vd =  0.9*BW;

    elseif GA <= 29
        BW =  1.44627e-06*(days^6)- 2.41190e-04*(days^5) + 1.58273e-02*(days^4)- ...
                  5.16579e-01*(days^3)+ 8.93645e00*(days^2)- 6.26752e01*(days)+ 1.11750e03;
        Ht =   3.61874e-08*(days^3)- 1.08198e-04*(days^2)+ 1.25024e-01*(days) + 3.53636e01;
        Vd =   0.85*BW;

    elseif GA <= 32
        BW =  -3.207719e-05*(days^6)+ 3.039153e-03*(days^5)- 1.080554e-01*(days^4) + ...
                  1.743312e00*(days^3)- 1.118843e01*(days^2) + 1.974714e01*(days)+ 1.382673e03;
        Ht =   3.727585e-08*(days^3) - 1.089800e-04*(days^2)+ 1.215347e-01*(days)+ 3.891467e01;
        Vd =   0.85*BW;

    elseif GA <= 36
        BW =  2.465457e-03*(days^4)- 1.807624e-01*(days^3)+ 5.953984e00*(days^2)-...
                  5.491700e01*(days)+ 1.958661e03;
        Ht =  - 3.401428e-09*(days^4) + 3.318122e-06*(days^3)- 1.173090e-03*(days^2) + ...
                 2.320859e-01*(days)+ 4.316323e01;
        Vd =  0.8*BW;

    else 
        BW =  - 1.262722e-08*(days^4) + 4.254248e-05*(days^3)- 5.008009e-02*(days^2) +...
                  3.086214e01*(days)+3.580496e03;
        Ht =   1.777345e-08*(days^3)- 6.086589e-05*(days^2)+ 8.644990e-02*(days)+ 5.148881e01;
        Vd =   0.75*BW;

    end

    TBW = 0.72 * BW;  % Females: 72% of BW
end

BW_kg = BW / 1000;  % Convert BW from g to kg
BSA = 0.024265 .* (BW_kg .^ 0.5378) .* (Ht .^ 0.3964);  % Haycock formula (m^2)

% GA/PNA/TH-aware parameters
sv_base_28w   = 22.0;
sv_slope_perW = -0.10;     % mL/m^2 per GA week
sv_sd         = 1.8;
sv_minmax     = [15 28];

sv_early_PNA48 = +1.0;
sv_TH_term     = -1.0;
sv_inotrope    = +1.0;

base = sv_base_28w + sv_slope_perW*(row.GA-28);
if row.PNA <= 2, base = base + sv_early_PNA48; end
if row.cooling==1 && row.GA>=36, base = base + sv_TH_term; end
if row.Ino==1, base = base + sv_inotrope; end
SV_index = normrnd(base, sv_sd);

% --- PDA effect (week-1), severity/PNA/BW dependent, with caps for VLBW
sev_gain = [0.00, 0.05, 0.10, 0.18];   % none/mild/moderate/severe (smaller)

if row.PDA==1 && row.PNA <= 6 && row.GA < 32
    % faster early decay
    decay = exp(-row.PNA/1.4);

    % BW-dependent damping: strongest <1.25 kg, fades by 2.0 kg
    bwkg = row.BW/1000;
    if bwkg < 1.0
        damp = 0.55;
    elseif bwkg < 1.5
        damp = 0.55 + 0.30*(bwkg - 1.0)/0.5;   % 0.55 -> 0.85 (1.0->1.5 kg)
    elseif bwkg < 2.0
        damp = 0.85 + 0.15*(bwkg - 1.5)/0.5;   % 0.85 -> 1.00 (1.5->2.0 kg)
    else
        damp = 1.00;
    end

    base_extra = sev_gain(row.PDA_sev+1) * decay * damp;

    % per-case hard cap to avoid tails (tighter for VLBW)
    max_extra = (bwkg < 1.25) * 0.12 + (bwkg >= 1.25) * 0.22;  % 0.12 vs 0.22
    extra = min(base_extra, max_extra);

    SV_index = SV_index * (1 + extra);
end


idx_vlbw = (row.BW/1000 < 1.25);
if any(idx_vlbw)
    mu_v = mean(SV_index(idx_vlbw),'omitnan');
    SV_index(idx_vlbw) = (SV_index(idx_vlbw)-mu_v)*(1.3/1.8) + mu_v;  % sd ~1.3 vs 1.8
end

% keep a conservative top cap
SV_index = min(max(SV_index, sv_minmax(1)), 30);

% Finish as before
row.SV_mL = SV_index .* BSA;                  % mL


HR = row.HR;

% ===== TH (cooled, term) — keep your Elstad-like time course =====
mu12 = 92;                 % target around 12 h
mu72 = 97;                 % target around 72 h
sdTH = 8;                  % spread during TH
ino_delta_TH = 15;         % inotrope effect during TH

% ===== non-TH baseline (GA-aware) =====
% Linear GA trend to match Hsu HR across 28–41 w:
% HR_mean ≈ a + b*GA with a≈203, b≈-1.92 ⇒ HR(28)≈149, HR(40)≈126
a_HR = 202.667; b_HR = -1.9167;
sdNT = 12;                 % spread for non-TH
ino_delta_NT = 5;          % small HR bump with inotropes when not cooled

if row.cooling==1 && row.GA>=36
    t = min(max(row.PNA_h,0),72);   % hours
    if t <= 12, mu = mu12;
    else,       mu = mu12 + (mu72-mu12)*(t-12)/(72-12);
    end
    if row.Ino==1, mu = mu + ino_delta_TH; end
    HR = round(normrnd(mu, sdTH));
    HR = min(max(HR, 80), 120);           % plausible TH bounds

else
    % Non-TH GA baseline (in bpm)
    mu = a_HR + b_HR * row.GA;

    % mild PNA effect (HR falls slightly after day 3; zero in 72–96h window)
    % Here we keep it minimal so Hsu window isn't distorted:
    if row.PNA > 5
        mu = mu - 3;   % tiny drop after day 5
    end

    % inotrope effect
    if row.Ino==1, mu = mu + ino_delta_NT; end

    HR = round(normrnd(mu, sdNT));
    HR = min(max(HR, 95), 170);           % non-TH bounds
end

row.HR = HR;
% --- 72–96 h non-TH HR pin to Hsu Table 2 means (additive offsets) -----
mask_hsu = (row.cooling==0) & (row.PNA>=3) & (row.PNA<5);

% GA bins as in Hsu (≤28, 29–30, 31–32, 33–34, 35–36, 37–38, 39–41)
ga_edges  = [0 28 30 32 34 36 38 41];
[g,~]     = discretize(row.GA, ga_edges);

% Hsu HR means by GA bin (bpm)
hsu_HR_mean = [149 145 142 142 136 131 126];

for k = 1:numel(hsu_HR_mean)
    idx = mask_hsu & (g==k);
    if any(idx)
        cur_mean = mean(row.HR,'omitnan');
        delta    = hsu_HR_mean(k) - cur_mean;   % additive correction
        row.HR = row.HR + round(delta);
    end
end

CO = (row.HR .* row.SV_mL) / 1000;  % Convert to L/min

% Update to table
row.PMA = PMA;
row.BW = BW;
row.Ht = Ht;
row.CO_L_min = CO;
row.TBW = TBW/1000;

% Vd amikacin
% da https://pmc.ncbi.nlm.nih.gov/articles/PMC5341529/#:~:text=They%20evaluated%20the%20pharmacokinetics%20in,was%205.3%20%C2%B1%203.3%20hours.

PNA=row.PNA;

    
    if GA <= 28                    % like PMA ≤29 when PNA is small
        if      PNA < 7,   Vd_per_kg = 0.414;   % their subgroup 1 mean
        elseif  PNA <= 28, Vd_per_kg = 0.472;   % subgroup 2 mean
        else               Vd_per_kg = 0.353;   % subgroup 3 mean (tiny n)
        end

    elseif GA <= 33                % like PMA 30–34
        if      PNA < 7,   Vd_per_kg = 0.462;   % subgroup 4 mean
        else               Vd_per_kg = 0.454;   % subgroup 5 mean
        end

    else                           % ga ≥ 34  (-> PMA ≥35, all PNA in paper)
        Vd_per_kg = 0.433;                    % subgroup 6 mean
    end

Vd_amikacin = Vd_per_kg .* BW_kg;
row.Vd_amikacin = Vd_amikacin;


updated_row = row;

end