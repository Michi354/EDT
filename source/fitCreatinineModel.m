function [bestFormula_char, lmeBest, meanRMSE, meanR2, figHandle] = fitCreatinineModel(selectedDataset, flag)

rng(42);
doTest = false;
% Determine output variable
colNames = selectedDataset.Properties.VariableNames;
idx = find(contains(colNames, 'rea', 'IgnoreCase', true));
if ~isempty(idx)
    selectedDataset.Properties.VariableNames{idx} = 'Creatinina';
end
outputVar_raw = 'Creatinina';
selectedDataset = selectedDataset(selectedDataset.(outputVar_raw) > 0, :); % avoid log(0)
selectedDataset.LogCreatinina = log(selectedDataset.(outputVar_raw));
outputVar = 'LogCreatinina';
idVar = 'PtID';

figHandle = []; 

% Remove rows with missing output
selectedDataset = selectedDataset(~isnan(selectedDataset.(outputVar)), :);
selectedDataset = selectedDataset(selectedDataset.BW~='x', :);
selectedDataset = selectedDataset(~isnan(selectedDataset.BW), :);
fprintf('Removed rows with NaN in %s. Remaining rows: %d\n', outputVar, height(selectedDataset));

% Covariate selection based on flag
switch flag
    case 0
        baselineCovs = {'BW', 'GA', 'PNA'};
        extraCovs = {'Ino', 'Ibu'};
        derivedTerms = {
            'PMA2', @(T) T.PMA.^2;
           
        };

        selectedDataset = nlterms(derivedTerms,selectedDataset);
        GA_bin = discretize(selectedDataset.GA, [24 40]);      % e.g., extreme, very, moderate, term
        BW_bin = discretize(selectedDataset.BW, [0 2500 5000]);
        groupingVariable = strcat("GA", string(GA_bin), "_BW", string(BW_bin));
        selectedDataset.(idVar) = string(selectedDataset.(idVar));
        
    case 1
        baselineCovs = {'BW',  'PNA','GA', 'aki'};
        extraCovs = {};
        derivedTerms = {
            'PMA2', @(T) T.PMA.^2;
            'BW2', @(T) T.BW.^2;
            'aki*PMA', @(T) double(T.aki).*T.PMA;
            'aki*BW', @(T) double(T.aki).*T.BW;
        };
        selectedDataset = nlterms(derivedTerms,selectedDataset);
        GA_bin = discretize(selectedDataset.GA, [24 26 27 29 30 32 33 36 37 42]);       % e.g., extreme, very, moderate, term
        BW_bin = discretize(selectedDataset.BW, [0 2500 5000]);
        groupingVariable = strcat("GA", string(GA_bin), "_BW", string(BW_bin));
        
        selectedDataset.PMA = selectedDataset.GA + round(selectedDataset.PNA/7);
        selectedDataset.(idVar) = string(selectedDataset.(idVar));

    case 2
        baselineCovs = {'BW', 'PNA','GA','cooling'};
        extraCovs = {'Ibu', 'Ino','asphyxia'};
        derivedTerms = {
            'PMA2', @(T) T.PMA.^2;
        };
        selectedDataset = nlterms(derivedTerms,selectedDataset);
        GA_bin = discretize(selectedDataset.GA, [24 26 27 29 30 32 33 36 37 42]);      % e.g., extreme, very, moderate, term
        BW_bin = discretize(selectedDataset.BW, [0 2500 5000]);
        groupingVariable = strcat("GA", string(GA_bin), "_BW", string(BW_bin));
        selectedDataset.(idVar) = string(selectedDataset.(idVar));

    case 3
        baselineCovs = {'BW', 'PNA','GA','cooling','aki'};
        extraCovs = {'Ibu', 'Ino'};
        derivedTerms = {
             'PMA2', @(T) T.PMA.^2;
            % 'BW2', @(T) T.BW.^2;
            % 'BWxAsphyxia', @(T) T.BW .* double(T.asphyxia);
        };
        selectedDataset = nlterms(derivedTerms,selectedDataset);
        GA_bin = discretize(selectedDataset.GA, [24 26 27 29 30 32 33 36 37 42]); 
        BW_bin = discretize(selectedDataset.BW, [0 1000 2000 3000 4000 5000]);
        groupingVariable = strcat("GA", string(GA_bin), "_BW", string(BW_bin));
        
        selectedDataset.PMA = selectedDataset.GA + round(selectedDataset.PNA/7);
        selectedDataset.(idVar) = string(selectedDataset.(idVar));


    otherwise
        error('Unknown flag value: %d', flag);
end

if ~ismember('PMA', baselineCovs)
    baselineCovs = ['PMA', baselineCovs];
end


% --- Identify valid covariates
allCovs = unique([baselineCovs, extraCovs, derivedTerms(:,1)']);
validCovs = {};
for i=1:numel(allCovs)
    if ismember(allCovs{i}, selectedDataset.Properties.VariableNames)
        v = selectedDataset.(allCovs{i});
        if isnumeric(v) 
            % && nanstd(double(v)) > 1e-8
            validCovs{end+1} = allCovs{i};
        elseif iscategorical(v) && numel(categories(removecats(v))) > 1
            validCovs{end+1} = allCovs{i};
        end
    end
end
fprintf('Valid covariates: %s\n', strjoin(validCovs, ', '));

if ~ismember('PMA', selectedDataset.Properties.VariableNames)
    selectedDataset.PMA = selectedDataset.GA + round(selectedDataset.PNA / 7);
end

% --- For flag=0: perform your existing K-fold CV ---
uniquePatients = unique(selectedDataset.(idVar));
kFolds = 5;
cv = cvpartition(groupingVariable, 'KFold', kFolds);

rmse_vals = NaN(kFolds,1);
mse_vals = NaN(kFolds,1);
R2_vals = NaN(kFolds,1);
R2_adj_vals = NaN(kFolds,1);
overallBestValR2 = -Inf;
overallBestFormula = '';
overallBestModel = [];
overallBest_T_val = [];

model_array = {};
model_formulas = {};
for fold=1:kFolds
    fprintf('\n-- Fold %d of %d --\n', fold, kFolds);
    trainIDs = selectedDataset.(idVar)(training(cv,fold));
    valIDs = selectedDataset.(idVar)(test(cv,fold));
    isTrain = ismember(selectedDataset.(idVar), trainIDs);
    isVal = ismember(selectedDataset.(idVar), valIDs);

    T_train = selectedDataset(isTrain,:);
    T_val = selectedDataset(isVal,:);

    fprintf('\n-- Training nID %d\tTesting nID %d\n', length(unique(string(trainIDs))), length(unique(string(valIDs))));

    if height(T_train)<2 || height(T_val)<2
        warning('Fold %d: Training or validation set too small.', fold);
        continue;
    end

    % Best random effects structure
    bestAIC = Inf; 
    %bestRE = sprintf('(1|%s)', idVar);
    bestRE = sprintf('(1 + PMA|%s)', idVar);
    if ismember('PNA', T_train.Properties.VariableNames)
        try
            testRE_PNA = '(1 + PNA|PtID)';
            testFormula_PNA = sprintf('%s ~ 1 + %s + %s', outputVar, strjoin(baselineCovs, ' + '), testRE_PNA);
            m_pna = fitlme(T_train, testFormula_PNA);
            if m_pna.ModelCriterion.AIC < bestAIC
                bestAIC = m_pna.ModelCriterion.AIC;
                bestRE = testRE_PNA;
            end
        catch
        end
    end
    for i=1:numel(baselineCovs)
        if ~ismember(baselineCovs{i}, T_train.Properties.VariableNames), continue; end
        testRE = sprintf('(1 + %s|%s)', baselineCovs{i}, idVar);
        testFormula = sprintf('%s ~ 1 + %s + %s', outputVar, strjoin(baselineCovs, ' + '), testRE);
        try
            m = fitlme(T_train, testFormula);
            if m.ModelCriterion.AIC < bestAIC
                bestAIC = m.ModelCriterion.AIC; 
                bestRE = testRE;
            end
        catch
        end
    end
    fprintf('Fold %d: Best random effects: %s\n', fold, bestRE);

    % Greedy forward selection
    currentCovs = baselineCovs;
    bestR2_train = -Inf;
    y_train = T_train.(outputVar);
    improvement = true;
    cnd = 0;

    while improvement
        improvement = false;
        candidates = setdiff(validCovs, currentCovs);
        bestNewR2 = bestR2_train;
        bestNewCov = '';
        for i=1:numel(candidates)
            testCovs = [currentCovs, candidates{i}];
            formulaStr = sprintf('%s ~ 1 + %s + %s', outputVar, strjoin(testCovs, ' + '), bestRE);
            try
            m = fitlme(T_train, formulaStr);
            cnd=cnd+1;
            model_array{fold, cnd} = m; model_formulas{fold, cnd} = formulaStr;
            y_pred = predict(m, T_train);
            toKeep = find(~isnan(y_pred));
            y_pred = y_pred(toKeep);
            y_train = y_train(toKeep);
            R2 = 1 - sum((y_train - y_pred).^2) / sum((y_train - mean(y_train)).^2);
            R2_adj = m.Rsquared.Adjusted;
            if R2_adj > bestNewR2
                bestNewR2 = R2_adj;
                bestNewCov = candidates{i};
                bestModel = m;
                bestFormula = formulaStr;
            end
            catch
            end
        end
        if ~isempty(bestNewCov) && bestNewR2 > bestR2_train
            currentCovs = [currentCovs, bestNewCov];
            bestR2_train = bestNewR2;
            improvement = true;
            fprintf('Added: %s (Train R²=%.6f)\n', bestNewCov, bestR2_train);
        end
    end

    % Validation performance
    try
        y_val = exp(T_val.(outputVar));  % convert log(Cr) back to Cr
        y_pred_val = exp(predict(bestModel, T_val)); % convert predicted log(Cr) back to Cr

        toKeep = find(~isnan(y_pred_val));
        y_pred_val = y_pred_val(toKeep);
        y_val = y_val(toKeep);
        rmse_vals(fold) = sqrt(mean((y_val - y_pred_val).^2));
        mse_vals(fold) = mean((y_val - y_pred_val).^2);
        R2_vals(fold) = 1 - sum((y_val - y_pred_val).^2) / sum((y_val - mean(y_val)).^2);
        R2_adj_vals(fold) = 1 - ((1 - R2_vals(fold)) * (height(T_val) - 1)) / (height(T_val) - numel(currentCovs) - 1);
        fprintf('Fold %d: RMSE=%.3f, R²=%.3f\n', fold, rmse_vals(fold), R2_vals(fold));
        if R2_adj_vals(fold) > overallBestValR2
            overallBestValR2 = R2_adj_vals(fold);
            overallBestFormula = bestFormula;
            overallBestModel = bestModel;
            overallBest_T_val = T_val;
        end
    catch
        warning('Fold %d: Validation prediction failed.', fold);
    end
end

% Final summary
bestFormula_char = overallBestFormula;
lmeBest = overallBestModel;
meanRMSE = mean(rmse_vals, 'omitnan');
meanMSE = mean(mse_vals, 'omitnan');
meanR2 = mean(R2_vals, 'omitnan');
meanR2_adj = mean(R2_adj_vals, 'omitnan');

fprintf('\n--- Final Summary ---\n');
fprintf('Best validation R²: %.3f\nFormula: %s\n', meanR2_adj, bestFormula_char);
fprintf('Best validation MSE: %.3f\n', meanMSE);
fprintf('Best validation RMSE: %.3f\n', meanRMSE);

if doTest
    for fold = 1:kFolds
        for i = 1:size(model_array,2)
            if not(isempty(model_array{fold,i}))
                if model_array{fold,i}.NumObservations == lmeBest.NumObservations
                [table,~] = compare(model_array{fold,i},lmeBest,nsim=50);
                fprintf('Model %s\n',model_formulas{fold,i})
                fprintf('pValue:  %.6f\n', table.pValue(2));
                end
            end
        end
    end
end

% Plot
if ~isempty(lmeBest) && exist('overallBest_T_val','var')
    y_true = overallBest_T_val.(outputVar);
    y_pred = predict(lmeBest, overallBest_T_val);
    figHandle = figure;
    scatter(y_true, y_pred, 'filled'); hold on;
    plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], 'k--', 'LineWidth', 1);
    xlabel('Observed'); ylabel('Predicted');
    title(sprintf('Validation: MSE=%.3f, R²=%.3f', meanMSE, meanR2_adj));
    grid on; hold off;
    mkdir './SavedFigures'
    saveas(gcf,fullfile('./SavedFigures',strcat('DiagnosticPlot_Model_',string(flag),'.png')))
end

end

function selectedDataset = nlterms(derivedTerms,selectedDataset)
    for i=1:size(derivedTerms,1)
        name = derivedTerms{i,1};
        fn = derivedTerms{i,2};
        try
            selectedDataset.(name) = fn(selectedDataset);
            selectedDataset = selectedDataset(~isnan(fn(selectedDataset)),:);
            % disp(size(selectedDataset,1))
        catch
        end
    end
end