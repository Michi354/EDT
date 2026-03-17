function co_affinity_quick_stats(T)
% Console-only “sanity” stats for CO/CI across common cohorts

need = {'CO_L_min','BW','BSA','GA','PNA','cooling'};
assert(all(ismember(need, T.Properties.VariableNames)), ...
    'co_affinity_quick_stats: missing columns: %s', ...
    strjoin(setdiff(need, T.Properties.VariableNames), ', '));

CO_Lmin = T.CO_L_min;
BW_kg   = T.BW/1000;
BSA     = T.BSA;

CO_mlkg = 1000*CO_Lmin./BW_kg;
CI      = CO_Lmin./BSA;

% 1) Doni-like
doni_mask = (T.GA>=37) & (T.PNA<=2) & (T.cooling==0);
fprintf('Doni-like (term, ≤48h, no TH): CO/kg median %0.0f; CI median %0.2f\n', ...
    median(CO_mlkg(doni_mask),'omitnan'), median(CI(doni_mask),'omitnan'));

% 2) Term under TH ≤48h
th_term48 = (T.GA>=37) & (T.PNA<=2) & (T.cooling==1);
fprintf('Term ≤48h under TH: CO/kg median %0.0f; CI median %0.2f\n', ...
    median(CO_mlkg(th_term48),'omitnan'), median(CI(th_term48),'omitnan'));

% 3) Hsu-like 72–96h, non-TH
hsu_mask = (T.PNA>=3) & (T.PNA<5) & (T.cooling==0);
fprintf('Hsu-like (72–96h, no TH): CI median %0.2f (target ~2.55)\n', ...
    median(CI(hsu_mask),'omitnan'));
end
