function elstad_hr_th_plots(T)
% Two boxcharts (~12h and ~72h) + inotrope slice, following Elstad-style windows

if ~ismember('PNA_h', T.Properties.VariableNames)
    assert(ismember('PNA', T.Properties.VariableNames), 'Missing PNA/PNA_h.');
    T.PNA_h = T.PNA * 24;
end

need = {'cooling','GA','HR','Ino','PNA_h'};
assert(all(ismember(need, T.Properties.VariableNames)), ...
    'elstad_hr_th_plots: missing columns: %s', ...
    strjoin(setdiff(need, T.Properties.VariableNames), ', '));

GA_min_TH  = 36;
ref12_med  = 92; ref12_CI = [88 96];
ref72_med  = 97; ref72_CI = [89 106];
winPNA     = [6 72];
win12      = [6 18];
win72      = [60 84];

isCoolRef = (T.cooling==1) & (T.GA>=GA_min_TH) & (T.PNA_h<=winPNA(2));
idx12 = isCoolRef & T.PNA_h>=win12(1) & T.PNA_h<=win12(2);
idx72 = isCoolRef & T.PNA_h>=win72(1) & T.PNA_h<=win72(2);

figure('Color','w'); tiledlayout(1,2,'Padding','compact');
nexttile;
boxchart(ones(sum(idx12),1), T.HR(idx12)); hold on;
yline(ref12_med,'--k','ref median 92');
patch(1+[-0.3 0.3 0.3 -0.3], [ref12_CI(1) ref12_CI(1) ref12_CI(2) ref12_CI(2)], ...
      [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.5);
ylabel('HR (bpm)'); title(sprintf('~12 h (n=%d)',sum(idx12)));
xlim([0.5 1.5]); grid on;

nexttile;
boxchart(ones(sum(idx72),1), T.HR(idx72)); hold on;
yline(ref72_med,'--k','ref median 97');
patch(1+[-0.3 0.3 0.3 -0.3], [ref72_CI(1) ref72_CI(1) ref72_CI(2) ref72_CI(2)], ...
      [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.5);
title(sprintf('~72 h (n=%d)',sum(idx72)));
xlim([0.5 1.5]); grid on;
sgtitle('Elstad-like windows with medians and 95% CI');

% Inotrope effect (12–24 h slice)
slice = (T.cooling==1) & (T.GA>=GA_min_TH) & T.PNA_h>=12 & T.PNA_h<=24;
groups = categorical(T.Ino(slice), [0 1], {'No inotrope','Inotrope'});
figure('Color','w'); boxchart(groups, T.HR(slice));
ylabel('HR (bpm)'); title('Inotrope effect on HR during TH (12–24 h)'); grid on;
end