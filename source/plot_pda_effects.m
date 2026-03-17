function plot_pda_effects(T, outPrefix)
% Creates the 3-panel PDA schematic and saves to PNG/PDF

if nargin<2 || isempty(outPrefix), outPrefix = 'pda_effect_schematic'; end

f = figure('Color','w','Position',[100 100 1200 420]);
t = tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');

% Panel A – prevalence with/without ibuprofen
PNA = 0:0.1:7; tau = 9; ibu_factor = 0.6;
ga_labels = {'\leq28 w','29–30 w','31–32 w','33–36 w','\geq37 w'};
p0        = [0.90,     0.80,      0.60,      0.35,      0.01];

prev_noIbu = zeros(numel(p0), numel(PNA));
prev_Ibu   = zeros(numel(p0), numel(PNA));
for j = 1:numel(p0)
    prev_noIbu(j,:) = p0(j) * exp(-PNA/tau);
    prev_Ibu(j,:)   = p0(j) * exp(-PNA/tau) * ibu_factor;
end

nexttile; hold on;
for j = 1:numel(p0), plot(PNA, prev_noIbu(j,:)*100, 'LineWidth', 2); end
for j = 1:numel(p0), plot(PNA, prev_Ibu(j,:)*100, '--', 'LineWidth', 1.5); end
grid on; xlabel('Postnatal age (days)'); ylabel('PDA prevalence (%)');
title('(A) PDA prevalence: GA bands, ibuprofen effect');
legend([strcat(ga_labels,' (no Ibu)'), strcat(ga_labels,' (+Ibu)')], 'Location','northeastoutside'); legend boxoff;

% Panel B – data-derived severity at Day1 vs Day5 from T.PDA_sev
nexttile;
assert(ismember('PDA', T.Properties.VariableNames) && ismember('PDA_sev', T.Properties.VariableNames), ...
       'plot_pda_effects: T must have PDA and PDA_sev');
isPDA  = T.PDA == 1;
d1mask = T.PNA >= 0 & T.PNA < 1 & isPDA;
d5mask = T.PNA >= 4 & T.PNA < 5 & isPDA;

sevCats = [1 2 3];
dist_D1 = arrayfun(@(s) mean(T.PDA_sev(d1mask) == s, 'omitnan'), sevCats);
dist_D5 = arrayfun(@(s) mean(T.PDA_sev(d5mask) == s, 'omitnan'), sevCats);
if ~any(d1mask), dist_D1(:) = NaN; end
if ~any(d5mask), dist_D5(:) = NaN; end

bar_data = 100 * [dist_D1; dist_D5];
bh = bar(categorical({'Day 1','Day 5'}), bar_data, 'stacked'); %#ok<NASGU>
ylabel('Distribution (%)'); title('(B) PDA severity D1 vs D5'); grid on; box on;
legend({'Mild','Moderate','Severe'},'Location','northeastoutside'); legend boxoff;

% Annotate
for ii = 1:2
    y0 = 0;
    for k = 1:3
        val = bar_data(ii,k);
        if isfinite(val) && val>3
            text(ii, y0 + val/2, sprintf('%.0f%%', val), 'HorizontalAlignment','center','FontSize',9);
        end
        y0 = y0 + (isfinite(val)*val);
    end
end

% Panel C – illustrative CO/kg with/without PDA
nexttile; hold on;
PNAc = 0:0.2:7;
COkg_vlbw_base = max(320 - 12*PNAc, 200);
COkg_term_base = max(190 - 3*PNAc, 150);
sev_gain = [0.00, 0.05, 0.10, 0.18]; sev_use = 2; decay = exp(-PNAc/1.4);
uplift_vlbw = 1 + sev_gain(sev_use+1) * decay * 0.70;
uplift_term = 1 + sev_gain(sev_use+1) * decay * 0.30;
plot(PNAc, COkg_vlbw_base, 'LineWidth', 2);
plot(PNAc, COkg_vlbw_base.*uplift_vlbw, '--', 'LineWidth', 2);
plot(PNAc, COkg_term_base, 'LineWidth', 2);
plot(PNAc, COkg_term_base.*uplift_term, '--', 'LineWidth', 2);
yline(200,'k:','Walther lower','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
yline(324,'k:','Walther upper','LabelHorizontalAlignment','left','LabelVerticalAlignment','top');
grid on; box on; xlabel('Postnatal age (days)'); ylabel('CO (mL/min/kg)');
title('(C) CO/kg with vs without PDA (illustrative)');
legend({'VLBW no PDA','VLBW + PDA','Term no PDA','Term + PDA'}, 'Location','northeast'); legend boxoff;

title(t,'PDA effects in the virtual cohort','FontWeight','bold');
drawnow;
exportgraphics(f, [outPrefix '.pdf'], 'ContentType','vector');
exportgraphics(f, [outPrefix '.png'], 'Resolution', 300);
fprintf('Saved: %s and %s\n', fullfile(pwd,[outPrefix '.pdf']), fullfile(pwd,[outPrefix '.png']));
end