function Tpred = prepare_data_baby(neonate, flag)

line   =neonate;

%% 1) Ensure grouping & outcome exist & are correct type

% % Make PtID categorical in both sets
% if ~iscategorical(line.(idVar)),   line.(idVar)   = categorical(line.(idVar));   end

%% 2) Augment both train & validation with derived covariates

if ismember('GA',line.Properties.VariableNames) && ismember('PNA',line.Properties.VariableNames) && ~ismember('PMA',line.Properties.VariableNames)
    line.PMA = line.GA + round(line.PNA/7);
end


if flag == 0   
    derived = {
        'PMA2', @(T) T.PMA.^2;
    };
elseif flag == 1   
    derived = {
            'PMA2', @(T) T.PMA.^2;
            'BW2', @(T) T.BW.^2;
            'aki*PMA', @(T) double(T.aki).*T.PMA;
            'aki*BW', @(T) double(T.aki).*T.BW;
        };
elseif flag == 2
     derived = {
        'PMA2', @(T) T.PMA.^2;
        'BW2', @(T) T.BW.^2;
        'PMAxBW', @(T) T.PMA .* T.BW;
        'PMAxCooling', @(T) T.PMA .* T.cooling;
        'BWxCooling', @(T) T.BW .* T.cooling;
        'BWxAsphyxia', @(T) T.BW .* double(T.asphyxia);
    };

elseif flag == 3
     derived = {
        'PMA2', @(T) T.PMA.^2;
        'BW2', @(T) T.BW.^2;
        'aki*PMA', @(T) double(T.aki).*T.PMA;
        'aki*BW', @(T) double(T.aki).*T.BW;
        'PMAxAsphyxia', @(T) T.PMA .* double(T.asphyxia);
        'PMAxCooling', @(T) T.PMA .* double(T.cooling);
    };
end


% Apply only those whose inputs exist
for i = 1:size(derived,1)
  name = derived{i,1}; fn = derived{i,2};
  % Check that fn(T) won't error by verifying all needed base names
  try
    line.(name) = fn(line);
  catch
    % skip if base vars missing
  end
end


%%

% GA_pop   = neonate.GA;
% PNA_pop  = neonate.PNA;
% PMA_pop = neonate.PMA;
% BW_pop   = neonate.BW;
% aki   = neonate.aki;
% Ibu = neonate.Ibu;
% Ino = neonate.Ino;
% asphyxia_pop   = neonate.cooling;

%% --- 2) Build the prediction table with *all* fixed‐effect vars ---
template = line(1, :);  
template.PtID = string(neonate.PtID);
% Now overwrite only the covariates from your new neonate:

% Initialize template as a struct

% BW and derived terms
if exist('BW_pop', 'var')
    template.BW = BW_pop;
    template.BW2 = BW_pop.^2;
end

% GA
if exist('GA_pop', 'var')
    template.GA = GA_pop;
end

% PMA and derived terms
if exist('PMA_pop', 'var')
    template.PMA = PMA_pop;
    template.PMA2 = PMA_pop.^2;
end

% PNA_Days
if exist('PNA_pop', 'var')
    template.PNA_Days = PNA_pop;
end

% aki (categorical or numeric)
if exist('aki', 'var')
    template.aki = aki;
end

% cooling
if exist('cooling', 'var')
    template.cooling = cooling;
end

% asphyxia
if exist('asphyxia_pop', 'var')
    template.asphyxia = asphyxia_pop;
end

% Ibu
if exist('Ibu', 'var')
    template.Ibu = Ibu;
end

% Ino
if exist('Ino', 'var')
    template.Ino = Ino;
end


% Gender (if needed)
if exist('Gender', 'var')
    template.Gender = Gender;
end

% Derived/interactions
if exist('PMA_pop', 'var') && exist('BW_pop', 'var')
    template.PMAxBW = PMA_pop .* BW_pop;
end
if exist('BW_pop', 'var') && exist('asphyxia_pop', 'var')
    template.BWxAsphyxia = BW_pop .* asphyxia_pop;
end
if exist('PMA_pop', 'var') && exist('asphyxia_pop', 'var')
    template.PMAxAsphyxia = PMA_pop .* asphyxia_pop;
end
if exist('PMA_pop', 'var') && exist('Ibu', 'var')
    template.PMAxIbu = PMA_pop .* Ibu;
end
if exist('GA_pop', 'var') && exist('cooling', 'var')
    template.GAxCool = GA_pop .* cooling;
end
if exist('GA_pop', 'var') && exist('Ino', 'var')
    template.GAxIno = GA_pop .* Ino;
end
if exist('PNA_pop', 'var') && exist('aki', 'var')
    template.PNA2xAki = (PNA_pop.^2) .* aki;
end
if exist('BW_pop', 'var') && exist('GA_pop', 'var')
    template.BWxGA = BW_pop .* GA_pop;
    template.BWxGA2 = (BW_pop .* GA_pop).^2;
end

Tpred = template;


end

