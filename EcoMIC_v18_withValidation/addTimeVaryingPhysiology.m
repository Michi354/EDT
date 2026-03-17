function [nEvents, tch] = addTimeVaryingPhysiology(m, sub, BW_kg, Vc)
% Add time-varying CL/flows (+ optional Kp_rem) as SimBiology events.
% Returns: nEvents (added) and tch (event times).

    % ---------- small helpers ----------
    hasParam = @(nm) ~isempty(sbioselect(m,'Type','parameter','Name',nm));
    fmt      = @(x) num2str(x,'%.15g');    % robust numeric formatting

    % ---------- sorted time grid ----------
    t = sub.TIME;  if ~isnumeric(t), t = str2double(string(t)); end
    t = double(t(:));
    [t, ord] = sort(t,'ascend');

    e  = []; c = []; GA = []; PNA = [];
    if ismember('eGFR',          sub.Properties.VariableNames), e   = sub.eGFR(ord);          end
    if ismember('CardiacOutput', sub.Properties.VariableNames), c   = sub.CardiacOutput(ord); end
    if ismember('GA',            sub.Properties.VariableNames), GA  = sub.GA(ord);            end
    if ismember('PNA',           sub.Properties.VariableNames), PNA = sub.PNA(ord);           end

    % convert to column doubles (if present)
    tocol = @(x) double(x(:));
    if ~isempty(e),   e   = tocol(e);   end
    if ~isempty(c),   c   = tocol(c);   end
    if ~isempty(GA),  GA  = tocol(GA);  end
    if ~isempty(PNA), PNA = tocol(PNA); end

    % ---------- forward-fill safely ----------
    if ~isempty(e)
        e = fillmissing(e,'previous');
        e = fillmissing(e,'nearest');
    end
    if ~isempty(c)
        c = fillmissing(c,'previous');
        c = fillmissing(c,'nearest');
    end
    if ~isempty(GA)
        GA = fillmissing(GA,'previous');
        GA = fillmissing(GA,'nearest');
    end
    if ~isempty(PNA)
        PNA = fillmissing(PNA,'previous');
        PNA = fillmissing(PNA,'nearest');
    end

    % ---------- change points (eGFR/CO/PNA) ----------
    change = true(size(t));
    if numel(t) > 1
        de = false(size(t)); dc = de; dp = de;
        if ~isempty(e),   de = [true; diff(e)   ~= 0]; end
        if ~isempty(c),   dc = [true; diff(c)   ~= 0]; end
        if ~isempty(PNA), dp = [true; diff(PNA) ~= 0]; end
        change = de | dc | dp;
    end

    tch = t(change);
    eg  = []; co = [];
    if ~isempty(e),   eg  = e(change);   end
    if ~isempty(c),   co  = c(change);   end
    if ~isempty(GA),  GA  = GA(change);  end
    if ~isempty(PNA), PNA = PNA(change); end

    % ---------- avoid events at obs/dose times ----------
    obsTimes  = [];
    doseTimes = [];
    if ismember('CONC', sub.Properties.VariableNames)
        obsTimes = double(sub.TIME(isfinite(sub.CONC)));
    end
    if ismember('DOSE', sub.Properties.VariableNames)
        doseTimes = double(sub.TIME(isfinite(sub.DOSE)));
    end
    anchors = unique([obsTimes(:); doseTimes(:)]);
    epsShift = 0.5; % h
    for j = 1:numel(tch)
        if tch(j) > 0 && any(abs(tch(j) - anchors) < 1e-9)
            tch(j) = max(0, tch(j) - epsShift);
        end
    end
    for j = 2:numel(tch)
        if tch(j) <= tch(j-1)
            tch(j) = tch(j-1) + 1e-3;
        end
    end

    % ---------- fixed volumes and baseline Kps ----------
    VH   = sbioselect(m,'Type','compartment','Name','Hearth');    VH = VH.Capacity;
    VK   = sbioselect(m,'Type','compartment','Name','Kidney');    VK = VK.Capacity;
    VR   = sbioselect(m,'Type','compartment','Name','Remainder'); VR = VR.Capacity;
    Kp_h = getParamOr(m,'PartitionCoefficientHearth',  0.25);
    Kp_k = getParamOr(m,'PartitionCoefficientKidney',  0.30);

    makeNonConstant(m,'VH');
    makeNonConstant(m,'VK');
    makeNonConstant(m,'VR');
    makeNonConstant(m,'Kp_h');
    makeNonConstant(m,'Kp_k');
    makeNonConstant(m,'CL');
    makeNonConstant(m,'BloodFlowRateHearth');
    makeNonConstant(m,'BloodFlowRateKidney');
    makeNonConstant(m,'BloodFlowRateRemainder');
    makeNonConstant(m,'PartitionCoefficientRemainder');



    % ---------- build events ----------
    nEvents = 0;
    for k = 1:numel(tch)
        % skip if key drivers are NaN
        if ~isempty(eg) && ~isfinite(eg(k)), continue; end
        if ~isempty(co) && ~isfinite(co(k)), continue; end

        exprs = {};

        % CL / kel from eGFR
        if ~isempty(eg)
            CL_L_h = eg(k) * BW_kg * 0.06;         % L/h
            kel_h  = CL_L_h / max(Vc, eps);
            if hasParam('CL')
                exprs{end+1} = ['CL='              fmt(CL_L_h)]; 
            end
            % if hasParam('EliminationRate')
            %     exprs{end+1} = ['EliminationRate=' fmt(kel_h) ]; 
            % end
        end

        % Flows from CO
        if ~isempty(co)
            CO_L_h = co(k) * 60;                   % L/h
            fh=0.05; fk=0.20; fr=1-fh-fk;
            s  = CO_L_h / max(CO_L_h*(fh+fk+fr), eps);
            Qh = CO_L_h*fh*s;  Qk = CO_L_h*fk*s;  Qr = CO_L_h*fr*s;
            if hasParam('BloodFlowRateHearth'),     exprs{end+1} = ['BloodFlowRateHearth='     fmt(Qh)]; end
            if hasParam('BloodFlowRateKidney'),     exprs{end+1} = ['BloodFlowRateKidney='     fmt(Qk)]; end
            if hasParam('BloodFlowRateRemainder'),  exprs{end+1} = ['BloodFlowRateRemainder='  fmt(Qr)]; end
        end

        % Optional Kp_rem from (GA,PNA)
        if ~isempty(GA) && ~isempty(PNA) && isfinite(GA(k)) && isfinite(PNA(k)) ...
                && hasParam('PartitionCoefficientRemainder')
            Vdpk = pickVdPerKg(GA(k), PNA(k));
            Vd   = Vdpk * BW_kg;
            num  = Vd - Vc - Kp_h*VH - Kp_k*VK;
            Kp_r = max(num / max(VR,eps), 0.05);
            exprs{end+1} = ['PartitionCoefficientRemainder=' fmt(Kp_r)];
        end

        if isempty(exprs), continue; end

        addevent(m, ['time >= ' fmt(tch(k))], exprs);  % cell array of assignments
        nEvents = nEvents + 1;
    end
end

% ---------- helpers ----------
function name = pickParam(m, candidates)
    name = "";
    for nm = string(candidates)
        p = sbioselect(m,'Type','parameter','Name',char(nm));
        if ~isempty(p), name = nm; return; end
    end
end

function cap = pickCap(m, candidates)
    cap = NaN;
    for nm = string(candidates)
        c = sbioselect(m,'Type','compartment','Name',char(nm));
        if ~isempty(c), cap = c.Capacity; return; end
    end
end

function makeNonConstant(m, names)
    for c = string(names)
        if strlength(c)==0, continue; end
        p = sbioselect(m,'Type','parameter','Name',char(c));
        if ~isempty(p)
            try, p.Constant = false; end %#ok<TRYNC>
        end
    end
end

function Vdpk = pickVdPerKg(GA, PNA)
    if GA <= 28
        if      PNA < 7,   Vdpk = 0.414;
        elseif  PNA <= 28, Vdpk = 0.472;
        else,              Vdpk = 0.353;
        end
    elseif GA <= 33
        if      PNA < 7,   Vdpk = 0.462;
        else,              Vdpk = 0.454;
        end
    else
        Vdpk = 0.433;
    end
end

function v = getParamOr(m, name, default)
    p = sbioselect(m,'Type','parameter','Name',name);
    if isempty(p), v = default; else, v = p.Value; end
end