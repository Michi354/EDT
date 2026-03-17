function [dataset, flag] = selectPipelineFromGUI()
    % Create a figure window
    f = uifigure('Name', 'Select Baby Condition', 'Position', [500, 500, 300, 300]);

    % Create a button group for the categories
    bg = uibuttongroup(f, 'Position', [50, 80, 200, 180], 'Title', 'Category');

    % Add radio buttons
    uiradiobutton(bg, 'Text', 'No Disease',          'Position', [10, 130, 150, 20]);
    uiradiobutton(bg, 'Text', 'Acute Kidney Injury', 'Position', [10, 100, 150, 20]);
    uiradiobutton(bg, 'Text', 'Asphyxia + TH',       'Position', [10, 70, 150, 20]);
    uiradiobutton(bg, 'Text', 'AKI + TH',            'Position', [10, 40, 150, 20]);


    % Create a button to confirm selection
    btn = uibutton(f, 'push', 'Text', 'Confirm', 'Position', [100, 40, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) confirmAndLoad(bg, f));

    % Wait for user to finish
    uiwait(f);

    % Load the dataset and flag saved in appdata during selection
    dataset = getappdata(0, 'dataset');
    flag = getappdata(0, 'flag');

    % Cleanup
    rmappdata(0, 'dataset');
    rmappdata(0, 'flag');
end

function confirmAndLoad(bg, f)
    % Retrieve selection
    selectedCategory = bg.SelectedObject.Text;
    fprintf('You selected: %s\n', selectedCategory);

    % File paths
    fn_amicool  = '.\dataverse_files\data as published amicool.xlsx';
    fn_creacool = '.\dataverse_files\dataset as published creacool.xlsx';
    fn_nonmem   = '.\dataverse_files\NONMEM.xlsx';

    % Initialize dataset and flag
    selectedDataset = [];
    flag = -1;

    % Load/merge based on category
    switch selectedCategory
        case 'No Disease'
            selectedDataset = readtable(fn_nonmem); 
            flag = 0;

        case 'Acute Kidney Injury'
            flag = 1;
            selectedDataset = readtable(fn_creacool);            

        case 'Asphyxia + TH'
            selectedDataset = readtable(fn_amicool);
            flag = 2;

        case 'AKI + TH'
            selectedDataset = enrichForDisease(fn_creacool, fn_amicool); %to have reasonable ensable values of patients with both aki and TH
            %selectedDataset = augmentDataset(selectedDataset,8);
            flag = 3;
        otherwise
            warning('Unknown selection! Returning empty dataset.');
    end
    

    % Save dataset and flag to appdata for retrieval
    setappdata(0, 'dataset', selectedDataset);
    setappdata(0, 'flag', flag);


    % Close the GUI window
    close(f);
end
