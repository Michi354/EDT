function guards = classical_schedules(row)
%CLASSICAL_SCHEDULES_FROM karel et al. (Table 1 report)

BW = row.BW;    % grams
PNA = row.PNA;  % days

guards = {};    

if PNA <= 14
    if BW <= 800
        guards{end+1} = struct('dose_mgkg',16,'tau_h',48,'inf_min',60);
    elseif (BW >800) && (BW <= 1200)
        guards{end+1} = struct('dose_mgkg',16,'tau_h',42,'inf_min',60);
    elseif (BW >1200) && (BW <= 2000)
        guards{end+1} = struct('dose_mgkg',15,'tau_h',36,'inf_min',60);
    elseif (BW >2000) && (BW <= 2800)
        guards{end+1} = struct('dose_mgkg',15,'tau_h',30,'inf_min',60);
    else
        guards{end+1} = struct('dose_mgkg',15,'tau_h',24,'inf_min',60);
    end

else 
    if BW <= 800
        guards{end+1} = struct('dose_mgkg',20,'tau_h',42,'inf_min',60);
    elseif (BW >800) && (BW <= 1200)
        guards{end+1} = struct('dose_mgkg',20,'tau_h',36,'inf_min',60);
    elseif (BW >1200) && (BW <= 2000)
        guards{end+1} = struct('dose_mgkg',18,'tau_h',30,'inf_min',60);
    elseif (BW >2000) && (BW <= 2800)
        guards{end+1} = struct('dose_mgkg',18,'tau_h',24,'inf_min',60);
    else
        guards{end+1} = struct('dose_mgkg',18,'tau_h',20,'inf_min',60);
    end
end

end
