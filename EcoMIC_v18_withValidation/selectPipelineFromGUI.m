function [dataset, flag] = selectPipelineFromGUI()
    f = uifigure('Name', 'Select Baby Condition', 'Position', [500, 500, 300, 300]);

    bg = uibuttongroup(f, 'Position', [50, 80, 200, 180], 'Title', 'Category');

    uiradiobutton(bg, 'Text', 'No Disease',          'Position', [10, 130, 150, 20]);
    uiradiobutton(bg, 'Text', 'Acute Kidney Injury', 'Position', [10, 100, 150, 20]);
    uiradiobutton(bg, 'Text', 'Asphyxia + TH',       'Position', [10, 70, 150, 20]);
    uiradiobutton(bg, 'Text', 'AKI + TH',            'Position', [10, 40, 150, 20]);


    btn = uibutton(f, 'push', 'Text', 'Confirm', 'Position', [100, 40, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) confirmAndLoad(bg, f));

    uiwait(f);

    dataset = getappdata(0, 'dataset');
    flag = getappdata(0, 'flag');

    rmappdata(0, 'dataset');
    rmappdata(0, 'flag');
end

function confirmAndLoad(bg, f)
    selectedCategory = bg.SelectedObject.Text;
    fprintf('You selected: %s\n', selectedCategory);

    fn_amicool  = '.\dataverse_files\data as published amicool.xlsx';
    fn_creacool = '.\dataverse_files\dataset as published creacool.xlsx';
    fn_nonmem   = '.\dataverse_files\NONMEM.xlsx';

    selectedDataset = [];
    flag = -1;

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
            selectedDataset = enrichForDisease(fn_creacool, fn_amicool);
            flag = 3;
        otherwise
            warning('Unknown selection! Returning empty dataset.');
    end
    

    setappdata(0, 'dataset', selectedDataset);
    setappdata(0, 'flag', flag);


    close(f);
end
