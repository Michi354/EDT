function simResults = createVirtualPop()
n_males = 500;
n_females = 500;
N = n_males + n_females;

bins_GA=[24 26 27 29 30 32 33 36 37 42];
GA = [  (26-24)*lhsdesign(n_males*0.1,1) + 24; 
        (29-27)*lhsdesign(n_males*0.1,1) + 27; 
        (32-30)*lhsdesign(n_males*0.2,1) + 30; 
        (36-33)*lhsdesign(n_males*0.3,1) + 33;  
        (42-37)*lhsdesign(n_males*0.3,1) + 37; 

        (26-24)*lhsdesign(n_females*0.1,1) + 24; 
        (29-27)*lhsdesign(n_females*0.1,1) + 27; 
        (32-30)*lhsdesign(n_females*0.2,1) + 30; 
        (36-33)*lhsdesign(n_females*0.3,1) + 33; 
        (42-37)*lhsdesign(n_females*0.3,1) + 37; 
     ]; % weeks

%PNA = 1+round(10 * lhsdesign(N,1).^4);
% Generate PNA in days (continuous, not rounded)
PNA_male = 28 * lhsdesign(n_males,1).^2;  % 0–28 d
PNA_fem  = 28 * lhsdesign(n_females,1).^2;
PNA      = [PNA_male; PNA_fem];

% PMA in weeks (keep as before, based on days)
PMA = GA + round(PNA/7);


BW = zeros(N,1);
Ht = zeros(N,1);
Vd = zeros(N,1);
HR = zeros(N,1);
% 
% simResults.ID = 1:N;
% simResults.neonates = zeros(N, 7); % GA, PNA, PMA, gender, TBW, Ht, Vd
% simResults.neonates(:,1) = round(GA);
% simResults.neonates(:,2) = round(PNA);
% simResults.neonates(:,3) = PMA;
% simResults.neonates(1:500,4) = 0;
% simResults.neonates(501:1000,4) = 1;

Gender = zeros(N,1);
Gender(1:500) = 0;
Gender(501:1000) = 1;


for n = 1:n_males
    days = PNA(n);
    if GA(n) <= 26
        BW(n) = 1.650186e-07*(days^6)- 3.972552e-05*(days^5) + 3.685725e-03*(days^4) ...
                 - 1.667271e-01*(days^3) + 4.012047e00*(days^2) - 3.705484e01*(days)+ 8.887167e02;
        Ht(n) =  8.905395e-04*(days^2)+ 4.139884e-02*(days)+ 3.375000e01;
        Vd(n) =  0.9*BW(n);

    elseif GA(n) <= 29
        BW(n) = 2.59724e-06*(days^6)- 4.23233e-04*(days^5) + 2.69872e-02*(days^4)...
                 - 8.49142e-01*(days^3)+ 1.38614e01*(days^2)- 9.52883e01*(days) + 1.18552e03;
        Ht(n) =  3.14905e-08*(days^3)- 1.05086e-04*(days^2)+ 1.29879e-01*(days)+ 3.52325e01;
        Vd(n) =  0.85*BW(n);

    elseif GA(n) <= 32
        BW(n) =  1.374450e-05*(days^6)- 1.525030e-03*(days^5)+ 6.816827e-02*(days^4)-...
                  1.566683e00*(days^3)+ 1.979696e01*(days^2)- 1.113418e02*(days)+ 1.565761e03;
        Ht(n) =   3.662092e-08*(days^3)- 1.099299e-04*(days^2)+1.254900e-01*(days)+ 3.901344e01;
        Vd(n) =   0.85*BW(n);

    elseif GA(n) <= 36
        BW(n) =  3.097386e-03*(days^4)- 1.981930e-01*(days^3)+ 5.740936e00*(days^2)-...
                  4.577231e01*(days)+ 1.978709e03;
        Ht(n) =   -3.501005e-09*(days^4)+ 3.405542e-06*(days^3)- 1.218821e-03*(days^2)+ 2.462250e-01*(days)+ 4.342661e01;
        Vd(n) =   0.8*BW(n);

    else 
        BW(n) =  1.076382e-11*(days^5)- 5.095864e-08*(days^4) + 9.191582e-05*(days^3)-...
                  7.825887e-02*(days^2)+ 3.746039e01*(days) + 3.579470e03;
        Ht(n) =   2.873236e-08*(days^3)- 8.491953e-05*(days^2) + 9.965299e-02*(days)+ 5.173054e01;
        Vd(n) =   0.75*BW(n);
    end
end

for n = 501:500+n_females
    days = PNA(n);
    if GA(n) <= 26
        BW(n) = 1.002907e-07*(days^6)- 2.415770e-05*(days^5)+ 2.234812e-03*(days^4)- 1.010066e-01*(days^3)+...
                 2.532719e00*(days^2) - 2.215122e01*(days) + 8.411017e02;
        Ht(n) =  7.419743e-04*(days^2)+ 5.619694e-02*(days) + 3.343850e01;
        Vd(n) =  0.9*BW(n);

    elseif GA(n) <= 29
        BW(n) =  1.44627e-06*(days^6)- 2.41190e-04*(days^5) + 1.58273e-02*(days^4)- ...
                  5.16579e-01*(days^3)+ 8.93645e00*(days^2)- 6.26752e01*(days)+ 1.11750e03;
        Ht(n) =   3.61874e-08*(days^3)- 1.08198e-04*(days^2)+ 1.25024e-01*(days) + 3.53636e01;
        Vd(n) =   0.85*BW(n);

    elseif GA(n) <= 32
        BW(n) =  -3.207719e-05*(days^6)+ 3.039153e-03*(days^5)- 1.080554e-01*(days^4) + ...
                  1.743312e00*(days^3)- 1.118843e01*(days^2) + 1.974714e01*(days)+ 1.382673e03;
        Ht(n) =   3.727585e-08*(days^3) - 1.089800e-04*(days^2)+ 1.215347e-01*(days)+ 3.891467e01;
        Vd(n) =   0.85*BW(n);

    elseif GA(n) <= 36
        BW(n) =  2.465457e-03*(days^4)- 1.807624e-01*(days^3)+ 5.953984e00*(days^2)-...
                  5.491700e01*(days)+ 1.958661e03;
        Ht(n) =  - 3.401428e-09*(days^4) + 3.318122e-06*(days^3)- 1.173090e-03*(days^2) + ...
                 2.320859e-01*(days)+ 4.316323e01;
        Vd(n) =  0.8*BW(n);

    else 
        BW(n) =  - 1.262722e-08*(days^4) + 4.254248e-05*(days^3)- 5.008009e-02*(days^2) +...
                  3.086214e01*(days)+3.580496e03;
        Ht(n) =   1.777345e-08*(days^3)- 6.086589e-05*(days^2)+ 8.644990e-02*(days)+ 5.148881e01;
        Vd(n) =   0.75*BW(n);

    end
end

% Build table, keep PNA in days
neonates_tbl = table(string(1:N)', round(GA), round(PNA), PMA, Gender, BW, Ht, Vd, ...
    'VariableNames', {'PtID','GA','PNA','PMA','Gender','BW','Ht','Vd'});

% Add PNA in hours (continuous, not rounded)
neonates_tbl.PNA_h = PNA * 24;
% simResults.neonates(:,5) = BW;
% simResults.neonates(:,6) = Ht;
% simResults.neonates(:,7) = Vd;

% AKI flag 
aki = zeros(N,1);
aki(1:n_males) = double(rand(N/2,1) < 0.0019);
aki(n_males+1:N) = double(rand(N/2,1) < 0.0135);
neonates_tbl.aki = aki;

% for i=1:500
%     simResults.neonates(:,8) = double(rand(N,1) < 0.0019); % nelle bambine è più probabile che nei bambini (7:1) con incidenza gen. 1.54%
% end 
% for i=501:1000
%     simResults.neonates(:,8) = double(rand(N,1) < 0.0135); 
% end

% COOLING flag ----> colonna 13

% Add COOLING flag
cooling = zeros(N,1);
for i=1:N      % Il cooling inizia <=6 ore e può durare fino a 72 ore (quindi al più il bambino sampled può essere safe con cooling a 1 fino a 72, in teoria 72+6 ma mettiamo così per sicurezza, cooling a 6 ore cominica ad essere rischioso)
    if neonates_tbl.GA(i) >= 36 && neonates_tbl.PNA_h(i) <= 72 && neonates_tbl.aki(i)==0
        cooling(i) = double(rand < 0.9);
    elseif neonates_tbl.GA(i) >= 36 && neonates_tbl.PNA_h(i) <= 72 && neonates_tbl.aki(i)==1
       cooling(i) = 1;
    else
       cooling(i) = 0;
    end

end

% For neonates with AKI = 1, assign cooling = 1 with 75% probability
% idxAki = (aki == 1);
% cooling(idxAki) = double(rand(sum(idxAki),1) < 0.75);

neonates_tbl.cooling = cooling;

% simResults.neonates(:,13) = zeros(N,1);
% for i=1:N
%     if simResults.neonates(i,1) >= 36 && PNA(i)*24 <=6
%         simResults.neonates(i,13) = double(rand < 0.9);  
%     else
%         simResults.neonates(i,13) = double(rand < 0.1);
%     end
% end

neonates_tbl.Ibu = double(rand(N,1) < 0.08);
neonates_tbl.Ino = double(rand(N,1) < 0.2);

%% --- PDA flag & severity (week-1, GA- and PNA-dependent) ---
pda = zeros(N,1);
pda_sev = zeros(N,1); % 0 none, 1 mild, 2 moderate, 3 severe

for i = 1:N
    ga  = neonates_tbl.GA(i);
    pna = neonates_tbl.PNA(i);  % days

    if ga < 28
        p0 = 0.9;
    elseif ga < 30
        p0 = 0.8;
    elseif ga < 32
        p0 = 0.6;
    elseif ga < 34
        p0 = 0.5;
    elseif ga < 36
        p0 = 0.35;
    elseif ga < 37
        p0 = 0.20;
    else
        p0 = 0.01; % term: uncommon after first 1–2 days
    end

    % Decay with PNA over the first week (spontaneous closure)
    p = p0 * exp(-pna/9);    % faster drop over week 1

    % Ibuprofen promotes closure → lower probability
    if neonates_tbl.Ibu(i)==1
        p = p * 0.6;
    end

    pda(i) = rand < p;

    if pda(i)
        % Severity distribution at presentation
        r = rand;
        if     r < 0.50, sev = 1;      % mild
        elseif r < 0.85, sev = 2;      % moderate
        else               sev = 3;    % severe
        end
        % After day 3, many ducts constrict → soften severity by one level
        if pna > 3, sev = max(1, sev-1); end

        % Ibuprofen can blunt severity (if not closing entirely)
        if neonates_tbl.Ibu(i)==1 && sev>1
            sev = sev - 1;
        end
        pda_sev(i) = sev;
    end
end

neonates_tbl.PDA     = pda;
neonates_tbl.PDA_sev = pda_sev;   % 0 none, 1 mild, 2 moderate, 3 severe


% --- Heart rate model (TH unchanged; non-TH GA-aware & Hsu-matched) ---
HR = zeros(N,1);

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

for i = 1:N
    if neonates_tbl.cooling(i)==1 && neonates_tbl.GA(i)>=36
        t = min(max(neonates_tbl.PNA_h(i),0),72);   % hours
        if t <= 12, mu = mu12;
        else,       mu = mu12 + (mu72-mu12)*(t-12)/(72-12);
        end
        if neonates_tbl.Ino(i)==1, mu = mu + ino_delta_TH; end
        HR(i) = round(normrnd(mu, sdTH));
        HR(i) = min(max(HR(i), 80), 120);           % plausible TH bounds

    else
        % Non-TH GA baseline (in bpm)
        mu = a_HR + b_HR * neonates_tbl.GA(i);

        % mild PNA effect (HR falls slightly after day 3; zero in 72–96h window)
        % Here we keep it minimal so Hsu window isn't distorted:
        if neonates_tbl.PNA(i) > 5
            mu = mu - 3;   % tiny drop after day 5
        end

        % inotrope effect
        if neonates_tbl.Ino(i)==1, mu = mu + ino_delta_NT; end

        HR(i) = round(normrnd(mu, sdNT));
        HR(i) = min(max(HR(i), 95), 170);           % non-TH bounds
    end
end

neonates_tbl.HR = HR;
% --- 72–96 h non-TH HR pin to Hsu Table 2 means (additive offsets) -----
mask_hsu = (neonates_tbl.cooling==0) & (neonates_tbl.PNA>=3) & (neonates_tbl.PNA<5);

% GA bins as in Hsu (≤28, 29–30, 31–32, 33–34, 35–36, 37–38, 39–41)
ga_edges  = [0 28 30 32 34 36 38 41];
[g,~]     = discretize(neonates_tbl.GA, ga_edges);

% Hsu HR means by GA bin (bpm)
hsu_HR_mean = [149 145 142 142 136 131 126];

for k = 1:numel(hsu_HR_mean)
    idx = mask_hsu & (g==k);
    if any(idx)
        cur_mean = mean(neonates_tbl.HR(idx),'omitnan');
        delta    = hsu_HR_mean(k) - cur_mean;   % additive correction
        neonates_tbl.HR(idx) = neonates_tbl.HR(idx) + round(delta);
    end
end


SV = zeros(N,1);        % mL

BW_kg = BW / 1000;  % Convert BW from g to kg
BSA = 0.024265 .* (BW_kg .^ 0.5378) .* (Ht .^ 0.3964);  % Haycock formula (m^2)

neonates_tbl.BSA = BSA;

%% --- Stroke volume / cardiac index physiology ---------------------------
% Needs: BW_kg, BSA, HR, GA, PNA, cooling, Ino already defined above.

%% --- Stroke volume / cardiac index physiology  ------
SV_index = zeros(N,1);  % mL/m^2

% GA/PNA/TH-aware parameters
sv_base_28w   = 22.0;
sv_slope_perW = -0.10;     % mL/m^2 per GA week
sv_sd         = 1.8;
sv_minmax     = [15 28];

sv_early_PNA48 = +1.0;
sv_TH_term     = -1.0;
sv_inotrope    = +1.0;

for i = 1:N
    base = sv_base_28w + sv_slope_perW*(neonates_tbl.GA(i)-28);
    if neonates_tbl.PNA(i) <= 2, base = base + sv_early_PNA48; end
    if neonates_tbl.cooling(i)==1 && neonates_tbl.GA(i)>=36, base = base + sv_TH_term; end
    if neonates_tbl.Ino(i)==1, base = base + sv_inotrope; end
    SV_index(i) = normrnd(base, sv_sd);
end

% --- PDA effect (week-1), severity/PNA/BW dependent, with caps for VLBW
sev_gain = [0.00, 0.05, 0.10, 0.18];   % none/mild/moderate/severe (smaller)

for i = 1:N
    if neonates_tbl.PDA(i)==1 && neonates_tbl.PNA(i) <= 6 && neonates_tbl.GA(i) < 32
        % faster early decay
        decay = exp(-neonates_tbl.PNA(i)/1.4);

        % BW-dependent damping: strongest <1.25 kg, fades by 2.0 kg
        bwkg = neonates_tbl.BW(i)/1000;
        if bwkg < 1.0
            damp = 0.55;
        elseif bwkg < 1.5
            damp = 0.55 + 0.30*(bwkg - 1.0)/0.5;   % 0.55 -> 0.85 (1.0->1.5 kg)
        elseif bwkg < 2.0
            damp = 0.85 + 0.15*(bwkg - 1.5)/0.5;   % 0.85 -> 1.00 (1.5->2.0 kg)
        else
            damp = 1.00;
        end

        base_extra = sev_gain(neonates_tbl.PDA_sev(i)+1) * decay * damp;

        % per-case hard cap to avoid tails (tighter for VLBW)
        max_extra = (bwkg < 1.25) * 0.12 + (bwkg >= 1.25) * 0.22;  % 0.12 vs 0.22
        extra = min(base_extra, max_extra);

        SV_index(i) = SV_index(i) * (1 + extra);
    end
end

% tighten dispersion for <1.25 kg to reduce SD
idx_vlbw = (neonates_tbl.BW/1000 < 1.25);
if any(idx_vlbw)
    mu_v = mean(SV_index(idx_vlbw),'omitnan');
    SV_index(idx_vlbw) = (SV_index(idx_vlbw)-mu_v)*(1.3/1.8) + mu_v;  % sd ~1.3 vs 1.8
end

% keep a conservative top cap
SV_index = min(max(SV_index, sv_minmax(1)), 30);

% Finish as before
SV = SV_index .* BSA;                  % mL
CO = (neonates_tbl.HR .* SV) / 1000;   % L/min
neonates_tbl.SV_mL     = SV;
neonates_tbl.CO_L_min  = CO;




%% === Pin CI to Hsu in 72–96 h (non-TH), GA-bin by GA-bin ===
mask_hsu = (neonates_tbl.cooling==0) & (neonates_tbl.PNA>=3) & (neonates_tbl.PNA<5);

% GA bins as in Hsu Table 2
ga_edges = [0 28 30 32 34 36 38 41];
[g,~]    = discretize(neonates_tbl.GA, ga_edges);

% Hsu CI means (L/min/m^2): ≤28, 29–30, 31–32, 33–34, 35–36, 37–38, 39–41
hsu_CI_mean = [2.31 2.45 2.69 2.54 2.64 2.49 2.60];

% Recover SV_index, compute current CI
SV_index = neonates_tbl.SV_mL ./ neonates_tbl.BSA;
CI_now   = neonates_tbl.CO_L_min ./ neonates_tbl.BSA;

for k = 1:numel(hsu_CI_mean)
    idx = mask_hsu & (g==k);
    if ~any(idx), continue; end
    cur = mean(CI_now(idx),'omitnan');    % match Hsu uses mean±SD
    if cur>0
        s = hsu_CI_mean(k)/cur;           % scale factor for this GA bin
        % Optional: keep scaling modest to avoid overfitting tiny bins
        s = max(min(s, 1.25), 0.75);
        SV_index(idx) = SV_index(idx) * s;
    end
end

% Recompute SV/CO after the pin
neonates_tbl.SV_mL    = SV_index .* neonates_tbl.BSA;
neonates_tbl.CO_L_min = (neonates_tbl.HR .* neonates_tbl.SV_mL) / 1000;


% Compute TBW in grams
TBW = zeros(N,1);
TBW(Gender == 0) = 0.79 * BW(Gender == 0);  % Males: 79% of BW
TBW(Gender == 1) = 0.72 * BW(Gender == 1);  % Females: 72% of BW

% Add to table (TBW in kg)
neonates_tbl.TBW = TBW/1000;

% Add Minimum Inhibitory Concentration (MIC)
% Ann Intensive Care. 2020 Jun 8;10:75. doi: 10.1186/s13613-020-00685-5
% Amikacin pharmacokinetic/pharmacodynamic in intensive care unit: a prospective database
micValue  = [2, 3, 4, 6, 8];  % tipici sono 2/4/6. 
micFreq   = [6, 2, 1, 0.5, 0.5];

% Sample MIC values from a discrete distribution using randsample
MIC = randsample(micValue , N, true , micFreq);
neonates_tbl.MIC = MIC';

%% Vd amikacin
% da https://pmc.ncbi.nlm.nih.gov/articles/PMC5341529/#:~:text=They%20evaluated%20the%20pharmacokinetics%20in,was%205.3%20%C2%B1%203.3%20hours.

% --- Amikacin Vd (L/kg) using GA × PNA bins, values from Hughes et al.
% JPPT 2017 (dopo suggerimento Karel)

Vd_per_kg = zeros(N,1);

for i = 1:N
    ga  = neonates_tbl.GA(i);   % weeks
    pna = neonates_tbl.PNA(i);  % days

    if ga <= 28                    % like PMA ≤29 when PNA is small
        if      pna < 7,   Vd_per_kg(i) = 0.414;   % their subgroup 1 mean
        elseif  pna <= 28, Vd_per_kg(i) = 0.472;   % subgroup 2 mean
        else               Vd_per_kg(i) = 0.353;   % subgroup 3 mean (tiny n)
        end

    elseif ga <= 33                % like PMA 30–34
        if      pna < 7,   Vd_per_kg(i) = 0.462;   % subgroup 4 mean
        else               Vd_per_kg(i) = 0.454;   % subgroup 5 mean
        end

    else                           % ga ≥ 34  (-> PMA ≥35, all PNA in paper)
        Vd_per_kg(i) = 0.433;                    % subgroup 6 mean
    end
end

BW_kg = neonates_tbl.BW/1000;
neonates_tbl.Vd_amikacin = Vd_per_kg .* BW_kg;

% Vd_amikacin = zeros(N,1);
% for i=1:N
%     if PMA(i)<=29 && PNA(i)<7
%         Vd_per_kg(i)=0.414;
%     elseif PMA(i)<=29 && PNA(i)<=28
%         Vd_per_kg(i)=0.472;
%     elseif PMA(i)<=29 && PNA(i)>28
%        Vd_per_kg(i)=0.353;  
%     elseif PMA(i)<=34 && PNA(i)<=7
%        Vd_per_kg(i)=0.462; 
%    elseif PMA(i)<=34 && PNA(i)>7
%        Vd_per_kg(i)=0.454;  
%     elseif PMA(i)>34
%        Vd_per_kg(i)=0.433; 
%     end
%     Vd_amikacin(i)=Vd_per_kg(i)*BW_kg(i);
% end
% neonates_tbl.Vd_amikacin = Vd_amikacin;

simResults.neonates = neonates_tbl;

figure
histogram(neonates_tbl.CO_L_min)
xlabel('Cardiac Output (L/min)')
ylabel('Number of neonates')
title('CO Distribution')
figure
% Define GA bins and group labels
edges = [24 28 32 36 40 44];
ga_bins = discretize(neonates_tbl.GA, edges);
labels = {'24–28', '28–32', '32–36', '36–40', '40–44'};

% Make boxplot
boxplot(neonates_tbl.CO_L_min, ga_bins)

% Set custom x-axis labels
set(gca, 'XTickLabel', labels)

xlabel('Gestational Age (weeks)')
ylabel('Cardiac Output (L/min)')
title('Cardiac Output by Gestational Age')

% % Step 1: Find neonates with Vd > 0.72
% idxHighVd = Vd_amikacin > 0.72;
% 
% % Step 2: Plot GA vs Vd for those patients
% figure;
% scatter(GA(idxHighVd), Vd_amikacin(idxHighVd), 50, 'r', 'filled');
% xlabel('Gestational Age (weeks)');
% ylabel('Vd Amikacin (L/kg)');
% title('Neonates with Vd > 0.72 L/kg');
% grid on;
% 
% % Step 3: Add horizontal line at 0.72 for reference
% yline(0.72, '--k', 'Vd = 0.72');
% 
% % Optional: Print these neonates' info
% disp('Neonates with Vd_amikacin > 0.72:');
% disp(simResults.neonates(idxHighVd, {'PtID','GA','PNA','BW','Vd_amikacin'}));


end


