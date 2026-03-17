function creacool_enriched_cleaned = enrichForDisease(fn_creacool, fn_amicool)

data_amicool_orig  = readtable(fn_amicool);
data_creacool_orig = readtable(fn_creacool);


outcomeVar='Creatinine_mg_dl_';
patientIDVar = 'PtID';
% --- Configuration ---
% Matching variables (must exist in both datasets after standardization)
matchingVars = {'GA', 'BW'};

% Variables to impute from amicool to creacool
% Ensure these are the correct names AFTER standardization in amicool
varsToImpute = {'Gender','cooling','asphyxia','Ibu','Ino','VEN','PMA'};

% Tolerances for nearest neighbor matching (if exact matches are too few)
tolerances.GA = 0;
tolerances.PNA_Days = 0;
tolerances.BW = 0;
% Set to 0 for exact match. Non-zero for range matching on numeric.

fprintf('--- Enriching Creacool Data by Matching with Amicool ---\n');

creacool_enriched = data_creacool_orig; 

actual_creacool_outcome_name = outcomeVar;
actual_creacool_id_name = 'nr';
if ismember(actual_creacool_outcome_name, creacool_enriched.Properties.VariableNames)
    creacool_enriched = renamevars(creacool_enriched, actual_creacool_outcome_name, outcomeVar);
end
if ismember(actual_creacool_id_name, creacool_enriched.Properties.VariableNames)
    creacool_enriched = renamevars(creacool_enriched, actual_creacool_id_name, patientIDVar);
end
if ismember(patientIDVar, creacool_enriched.Properties.VariableNames) && isnumeric(creacool_enriched.(patientIDVar))
    creacool_enriched.(patientIDVar) = categorical(creacool_enriched.(patientIDVar));
end
fprintf('Creacool data loaded. Standardized names: %s -> %s, %s -> %s.\n', ...
    actual_creacool_outcome_name, outcomeVar, actual_creacool_id_name, patientIDVar);

% --- 2. Load and Standardize amicool data ---
opts_amicool = detectImportOptions(fn_amicool);
% Dynamically set categorical types for varsToImpute if they are read as non-categorical
% Get expected types from a dummy read or prior knowledge
temp_opts_for_type_check = detectImportOptions(fn_amicool); % Read options again for safety
original_amicool_varnames = temp_opts_for_type_check.VariableNames;
amicool_vars_to_set_categorical = {};
for v_idx = 1:length(varsToImpute)
    var_to_check = varsToImpute{v_idx};
    % Find the original name of var_to_check in amicool file if it's different
    % For simplicity, assuming varsToImpute uses the names as they appear in amicool file
    % OR that they will be standardized correctly below.
    % This part is tricky if varsToImpute are already standardized names.
    % Let's assume varsToImpute are the names *in the amicool file* for type detection.
    if ismember(var_to_check, original_amicool_varnames)
        idx = find(strcmp(original_amicool_varnames, var_to_check));
        if ~strcmp(temp_opts_for_type_check.VariableTypes{idx}, 'categorical') && ...
           (any(strcmp(var_to_check, {'Gender', 'cooling', 'asphyxia', 'Ibu', 'Ino', 'VEN'}))) % Add other known categoricals
            amicool_vars_to_set_categorical{end+1} = var_to_check; %#ok<AGROW>
        end
    end
end
if ~isempty(amicool_vars_to_set_categorical)
    fprintf('Setting the following amicool variables to categorical upon import: %s\n', strjoin(amicool_vars_to_set_categorical, ', '));
    opts_amicool = setvartype(opts_amicool, amicool_vars_to_set_categorical, 'categorical');
end

data_amicool_standardized = data_amicool_orig; 

actual_amicool_outcome_name = 'Creatinine_mg_dl'; 
actual_amicool_id_name = 'PtID'; 
if ismember(actual_amicool_outcome_name, data_amicool_standardized.Properties.VariableNames) && ~strcmp(actual_amicool_outcome_name, outcomeVar)
    data_amicool_standardized = renamevars(data_amicool_standardized, actual_amicool_outcome_name, outcomeVar);
end
if ismember(actual_amicool_id_name, data_amicool_standardized.Properties.VariableNames) && ~strcmp(actual_amicool_id_name, patientIDVar)
    data_amicool_standardized = renamevars(data_amicool_standardized, actual_amicool_id_name, patientIDVar);
end
if ismember(patientIDVar, data_amicool_standardized.Properties.VariableNames) && isnumeric(data_amicool_standardized.(patientIDVar))
    data_amicool_standardized.(patientIDVar) = categorical(data_amicool_standardized.(patientIDVar));
end

missing_amicool_vars = setdiff([varsToImpute, matchingVars], data_amicool_standardized.Properties.VariableNames);
if ~isempty(missing_amicool_vars)
    error('Following variables are missing from amicool data after standardization: %s. Check `varsToImpute` and `matchingVars` names.', strjoin(missing_amicool_vars, ', '));
end
fprintf('Amicool data loaded. Standardized names if necessary.\n');

% --- 3. Prepare creacool_enriched table for new columns ---
for k = 1:numel(varsToImpute)
    varName = varsToImpute{k};
    if ~ismember(varName, data_amicool_standardized.Properties.VariableNames)
        error('Variable "%s" specified in varsToImpute does not exist in standardized amicool data.', varName);
    end
    
    source_col_in_amicool = data_amicool_standardized.(varName);

    if iscategorical(source_col_in_amicool)
        cats = categories(source_col_in_amicool);
        isOrd = isordinal(source_col_in_amicool);
        if isempty(cats) % if amicool column is all <undefined> or empty itself
             creacool_enriched.(varName) = categorical(repmat(missing, height(creacool_enriched), 1));
        else
             creacool_enriched.(varName) = categorical(repmat(missing, height(creacool_enriched), 1), cats, 'Ordinal', isOrd);
        end
    elseif isnumeric(source_col_in_amicool)
        creacool_enriched.(varName) = NaN(height(creacool_enriched), 1);
    elseif islogical(source_col_in_amicool)
        creacool_enriched.(varName) = false(height(creacool_enriched),1); % Default to false, or use NaN if preferred
    elseif iscellstr(source_col_in_amicool) || isstring(source_col_in_amicool)
        creacool_enriched.(varName) = repmat({''}, height(creacool_enriched), 1); % For cell strings
        if isstring(source_col_in_amicool)
            creacool_enriched.(varName) = strings(height(creacool_enriched),1); % For string array
        end
    else
        warning('Unsupported data type for variable %s in amicool. Column not added properly.', varName);
        % Or error out: error('Unsupported data type for variable %s', varName);
    end
end
fprintf('Added empty columns for imputed vars to creacool_enriched.\n');

% --- 4. Perform Matching and Imputation ---
unmatched_creacool_indices = [];
fprintf('Starting matching process for %d creacool measurements...\n', height(creacool_enriched));

for i = 1:height(creacool_enriched)
    if mod(i, 50) == 0 || i == height(creacool_enriched) % Print progress less frequently
        fprintf('Processing creacool measurement %d of %d...\n', i, height(creacool_enriched));
    end

    creacool_subject_match_values = struct();
    valid_creacool_match_data = true;
    for k_match=1:length(matchingVars)
        var = matchingVars{k_match};
        if ~ismember(var, creacool_enriched.Properties.VariableNames)
            warning('Matching variable %s not found in creacool_enriched. Skipping subject %d for this var.', var, i);
            valid_creacool_match_data = false; break;
        end
        val = creacool_enriched.(var)(i);
        if (isnumeric(val) && isnan(val)) || (iscategorical(val) && isundefined(val))
            warning('Matching variable %s is NaN/undefined for creacool subject %d. Cannot match this subject.', var, i);
            valid_creacool_match_data = false; break;
        end
        creacool_subject_match_values.(var) = val;
    end

    if ~valid_creacool_match_data
        unmatched_creacool_indices(end+1) = i;
        continue;
    end

    amicool_indices = true(height(data_amicool_standardized), 1); 

    for k_match=1:length(matchingVars)
        var = matchingVars{k_match};
        cc_val = creacool_subject_match_values.(var);
        am_vals = data_amicool_standardized.(var);
        tol = tolerances.(var);

        if tol == 0 
            if iscategorical(cc_val)
                amicool_indices = amicool_indices & (am_vals == cc_val);
            else 
                amicool_indices = amicool_indices & (abs(am_vals - cc_val) < 1e-6); % Check for numeric equality robustly
            end
        else 
             if iscategorical(cc_val) 
                 warning('Tolerance matching for categorical var %s not implemented, using exact.',var);
                 amicool_indices = amicool_indices & (am_vals == cc_val);
             else 
                amicool_indices = amicool_indices & (am_vals >= (cc_val - tol) & am_vals <= (cc_val + tol));
             end
        end
    end
    
    matched_amicool_subset = data_amicool_standardized(amicool_indices, :);

    if height(matched_amicool_subset) > 0
        for j = 1:numel(varsToImpute)
            varName = varsToImpute{j};
            source_column_data = matched_amicool_subset.(varName);

            if iscategorical(creacool_enriched.(varName))
                if all(isundefined(source_column_data))
                    creacool_enriched.(varName)(i) = undefined;
                else
                    defined_values_in_source = source_column_data(~isundefined(source_column_data));
                    if ~isempty(defined_values_in_source)
                        try
                            mode_val = mode(defined_values_in_source);
                            creacool_enriched.(varName)(i) = mode_val;
                        catch ME_mode
                            warning('Could not compute mode for categorical %s for creacool subject %d. Details: %s. Setting to undefined.', varName, i, ME_mode.message);
                            creacool_enriched.(varName)(i) = undefined;
                        end
                    else
                        creacool_enriched.(varName)(i) = undefined;
                    end
                end
            elseif isnumeric(creacool_enriched.(varName))
                if all(isnan(source_column_data))
                    creacool_enriched.(varName)(i) = NaN;
                else
                    creacool_enriched.(varName)(i) = mean(source_column_data, 'omitnan');
                end
            elseif islogical(creacool_enriched.(varName))
                if isempty(source_column_data) || all(source_column_data == false & isnan(double(source_column_data))) % check if all are effectively missing/false
                     creacool_enriched.(varName)(i) = false; % Default if no info or all false
                else
                    % mode() on logical might be tricky if all are same.
                    % A simple majority vote:
                    if sum(source_column_data) >= (numel(source_column_data) / 2)
                        creacool_enriched.(varName)(i) = true;
                    else
                        creacool_enriched.(varName)(i) = false;
                    end
                end
            elseif iscellstr(creacool_enriched.(varName)) || isstring(creacool_enriched.(varName)) % Handles cell strings and string arrays
                 if isempty(source_column_data) || (iscell(source_column_data) && all(cellfun('isempty', source_column_data))) || (isstring(source_column_data) && all(source_column_data == ""))
                    if isstring(creacool_enriched.(varName))
                        creacool_enriched.(varName)(i) = "";
                    else
                        creacool_enriched.(varName){i} = '';
                    end
                 else
                    % Ensure source_column_data is cellstr for unique/accumarray if it was string
                    if isstring(source_column_data)
                        source_column_data_cell = cellstr(source_column_data(source_column_data~="")); % convert non-empty strings to cell
                    else % already cellstr
                        source_column_data_cell = source_column_data(~cellfun('isempty', source_column_data));
                    end

                    if ~isempty(source_column_data_cell)
                        [unique_strings, ~, string_idx] = unique(source_column_data_cell);
                        counts = accumarray(string_idx, 1);
                        [~, max_idx_str] = max(counts);
                        if isstring(creacool_enriched.(varName))
                             creacool_enriched.(varName)(i) = string(unique_strings{max_idx_str});
                        else
                             creacool_enriched.(varName){i} = unique_strings{max_idx_str};
                        end
                    else
                        if isstring(creacool_enriched.(varName))
                            creacool_enriched.(varName)(i) = "";
                        else
                            creacool_enriched.(varName){i} = '';
                        end
                    end
                 end
            end
        end
    else
        unmatched_creacool_indices(end+1) = i;
    end
end

fprintf('Matching and imputation complete.\n');
if ~isempty(unmatched_creacool_indices)
    fprintf('%d out of %d creacool subjects could not be matched based on current criteria.\n', ...
            numel(unmatched_creacool_indices), height(creacool_enriched));
else
    fprintf('All %d creacool subjects were processed (either matched or imputation attempted).\n', height(creacool_enriched));
end
fprintf('--- End of Enrichment Process ---\n');

fprintf('--- Removing rows with missing values in imputed columns from creacool_enriched ---\n');
creacool_enriched_cleaned = creacool_enriched; % Work with a copy

rows_to_keep = true(height(creacool_enriched_cleaned), 1);

for k_impute_check = 1:numel(varsToImpute)
    varName_check = varsToImpute{k_impute_check};
    if ~ismember(varName_check, creacool_enriched_cleaned.Properties.VariableNames)
        warning('Variable %s (for NaN check) not found in creacool_enriched_cleaned. Skipping.', varName_check);
        continue;
    end
    
    current_column_to_check = creacool_enriched_cleaned.(varName_check);
    
    if iscategorical(current_column_to_check)
        rows_to_keep = rows_to_keep & ~isundefined(current_column_to_check);
    elseif isnumeric(current_column_to_check)
        rows_to_keep = rows_to_keep & ~isnan(current_column_to_check);
    elseif islogical(current_column_to_check)
        % Assuming logical false is valid, not missing. If NaNs were used for missing logicals:
        % rows_to_keep = rows_to_keep & ~isnan(double(current_column_to_check));
    elseif iscellstr(current_column_to_check)
        rows_to_keep = rows_to_keep & ~cellfun('isempty', current_column_to_check);
    elseif isstring(current_column_to_check)
        rows_to_keep = rows_to_keep & (current_column_to_check ~= "");
    end
end

num_rows_before_final_removal = height(creacool_enriched_cleaned);
creacool_enriched_cleaned = creacool_enriched_cleaned(rows_to_keep, :);
num_rows_after_final_removal = height(creacool_enriched_cleaned);
num_removed_due_to_missing_imputed = num_rows_before_final_removal - num_rows_after_final_removal;

fprintf('Number of rows in enriched table before final cleaning: %d\n', num_rows_before_final_removal);
fprintf('Number of rows after removing those with missing/undefined imputed values: %d\n', num_rows_after_final_removal);
fprintf('%d rows were removed due to missing values in the imputed columns.\n', num_removed_due_to_missing_imputed);

fprintf('--- End of Enrichment and Cleaning Process ---\n');
end