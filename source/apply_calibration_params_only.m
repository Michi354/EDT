function apply_calibration_params_only(modelObj, row, p, theta_hat, opts)
    if nargin<5, opts = struct; end
    if ~isfield(opts,'y0'), opts.y0 = [0;0;0;0;0;0; 4.1967e5; 125901]; end

    m  = modelObj;
    cs = getconfigset(m,'active');
    % You can leave StopTime/tolerances as-is; no sim here
    set(cs.SolverOptions, 'AbsoluteTolerance', 1e-9, 'RelativeTolerance', 1e-6);

    % ---------- Physiology ----------
    setCompCap(m,'Hearth',       p.Hearth,       'liter');
    setCompCap(m,'Kidney',       p.Kidney,       'liter');
    setCompCap(m,'Remainder',    p.Remainder,    'liter');
    setCompCap(m,'Serum_blood',  p.Serum_blood,  'liter');

    setP(m,'BloodFlowRateHearth',    p.BloodFlowRateHearth);
    setP(m,'BloodFlowRateKidney',    p.BloodFlowRateKidney);
    setP(m,'BloodFlowRateRemainder', p.BloodFlowRateRemainder);

    setP(m,'PartitionCoefficientHearth',    p.PartitionCoefficientHearth);
    setP(m,'PartitionCoefficientKidney',    p.PartitionCoefficientKidney);
    setP(m,'PartitionCoefficientRemainder', p.PartitionCoefficientRemainder);

    setP(m,'CL', p.CL);
    maybeSet(m,'kf',  getfielddef(opts,'kf',  p.kf, 1));
    maybeSet(m,'fu',  getfielddef(opts,'fu',  p.fu, 0.9));
    maybeSet(m,'aRS', getfielddef(p,'aRS', 0.5));
    maybeSet(m,'aSR', getfielddef(p,'aSR', 0.5));
    maybeSet(m,'K',   getfielddef(p,'K',   2.32e9));

    % ---------- PD theta ----------
    % [EC50s0 EC50r0 c rS Emax rR ku sigma d]
    setP(m,'EC50s0', theta_hat(1));
    setP(m,'EC50r0', theta_hat(2));
    setP(m,'c',      theta_hat(3));
    setP(m,'rS',     theta_hat(4));
    setP(m,'Emax',   theta_hat(5));
    setP(m,'rR',     theta_hat(6));
    setP(m,'ku',     theta_hat(7));
    setP(m,'sigma',  theta_hat(8));
    setP(m,'d',      theta_hat(9));

    % ---------- Initial conditions ----------
    y0 = opts.y0;
    setIC(m,'Cb',             y0(1));
    setIC(m,'Aminoglycoside', y0(2));
    setIC(m,'Ch',             y0(3));
    setIC(m,'Ck',             y0(4));
    setIC(m,'Cr',             y0(5));
    setIC(m,'u',              y0(6));
    setIC(m,'Ns',             y0(7));
    setIC(m,'Nr',             y0(8));
end

% ---------- helpers ----------
function setP(model, name, val)
    p = sbioselect(model,'Type','parameter','Name',name);
    if isempty(p), return; end
    p.Value = val;
end

function setIC(model, speciesName, val)
    s = sbioselect(model,'Type','species','Name',speciesName);
    if isempty(s), return; end
    s.InitialAmount = val;
end

function setCompCap(model, compName, cap, units)
    c = sbioselect(model,'Type','compartment','Name',compName);
    if isempty(c), return; end
    c.Capacity = cap;
    try, c.CapacityUnits = units; end
end

function maybeSet(model, parName, val)
    try, setP(model,parName,val); catch, end
end

function out = getfielddef(S, name, defaultVal, altDefault)
    if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
        out = S.(name);
    elseif nargin>=4
        out = altDefault;
    else
        out = defaultVal;
    end
end
