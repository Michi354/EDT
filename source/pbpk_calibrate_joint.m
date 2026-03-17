function [theta_hat, p] = pbpk_calibrate_joint(row, opts)
% Fit PD/eco parameters only to make classical schedules behave well.
% No paper targets are used.
%
% theta = [EC50s0, EC50r0, c, rS, Emax, rR, ku, sigma, d]
%
% opts.guard_schedules : {struct('dose_mgkg',..,'tau_h',..,'inf_min',60), ...}
% opts.duration_h      : simulation horizon (default 72h)
% opts.alpha           : potency conversion (gentamicin -> amikacin) for schedule conversion (default 4)
% opts.inf_min_ami     : infusion duration (min) (default 60)
% opts.units           : 'molecule_per_L' (default) or 'CFU_per_mL'

if nargin < 2, opts = struct; end
    if ~isfield(opts,'duration_h'), opts.duration_h=72; end
    if ~isfield(opts,'fu'),         opts.fu=1; end
    if ~isfield(opts,'kf'),         opts.kf=1; end
    if ~isfield(opts,'inf_min_ami'),opts.inf_min_ami=60; end
    if ~isfield(opts,'units'),      opts.units='molecule_per_L'; end
    if ~isfield(opts,'alpha'),      opts.alpha=4; end % <-- only for gentamicin->amikacin eq
    if ~isfield(opts,'guard_schedules')||isempty(opts.guard_schedules)
        opts.guard_schedules = classical_schedules(row);
    end

    % scale factor
    if strcmpi(opts.units,'molecule_per_L')
        opts.N_scale=1000;
    else
        opts.N_scale=1;
    end

    % free parameters (EC50r0 now free)
    th0 = [8, 8, 0.12, 0.9, 6, 0.70, 0.8, 0.18, 0.08];
    lb = [6, 6, 0.02, 0.4, 1.0, 0.10, 0.2, 0.05, 0.05];
    ub = [40,40, 0.60, 1.5,10.0, 0.70, 1.0, 0.30, 0.18];

    if isfield(opts,'theta0'), th0=opts.theta0; end
    if isfield(opts,'lb'),     lb=opts.lb; end
    if isfield(opts,'ub'),     ub=opts.ub; end

    x0 = min(max(th0,lb),ub);

    % guard-only residual
    BW_kg  = row.BW/1000;
    if ~isfinite(BW_kg) || BW_kg<=0
            if ismember('BW',row.Properties.VariableNames) && isfinite(row.BW(1))
                BW_kg = row.BW(1)/1000;
            elseif ismember('CW',row.Properties.VariableNames) && isfinite(row.CW(1))
                BW_kg = row.CW(1)/1000;
            else
                BW_kg = 3.63;
            end
    end
    row.BW_kg = BW_kg;
    p_base = build_phys_from_row(row,opts);
    
   

    function r = guard_resid(th)
        Ns0 = 4.1967e5 * opts.N_scale;
        Nr0 = 125901   * opts.N_scale;
        y0  = [0;0;0;0;0;0;Ns0;Nr0];
        N0  = Ns0 + Nr0;
        r=[];
        for g=1:numel(opts.guard_schedules)
            G=opts.guard_schedules{g};
            t_end=max(opts.duration_h,3*G.tau_h);
            times=0:G.tau_h:t_end;

            [t,y]=simulate_pbpk_one(p_base,BW_kg,G.dose_mgkg,times,t_end,y0,th,G.inf_min);
            Amino=y(:,2); Ns=y(:,7); Nr=y(:,8); N=Ns+Nr;

            % (1) Peak total ≤ baseline
            r=[r; max(0,log10(max(N)+1e-12)-log10(N0))];
            % (2) End ≤ baseline
            r=[r; max(0,log10(N(end)+1e-12)-log10(N0))];
            % (3) Resistant peak ≤ 10× initial Nr
            r=[r; max(0,log10(max(Nr)+1e-12)-log10(10*Nr0+1e-12))];
            % (4) Avoid full early Ns wipeout
            early=t<=min(12,G.tau_h);
            if any(early)
                Ns_floor=3e6; 
                minNs=max(min(Ns(early)),1e-12);
                r=[r; max(0,log10(Ns_floor/minNs))];
            end
            % (5) During on-drug, R should not be more suppressed than S
            dA=gradient(Amino,t);
            on=(Amino>0 & dA>=-1e-8);
            if any(on)
                dNs=gradient(Ns,t); dNr=gradient(Nr,t);
                adv=mean(dNr(on)-dNs(on));
                r=[r; max(0,-adv)];
            end
            % (6) Limit early drop in total burden (too-abrupt kill)
            early_win = t <= min(12, G.tau_h);
            if any(early_win)
                N_start = max(N(find(early_win,1,'first')), 1e-12);
                N_min   = max(min(N(early_win)), 1e-12);
                drop_log = log10(N_start) - log10(N_min);
                drop_cap = 0.5;                % allow up to ~0.8 log10 drop early
                r = [r; max(0, drop_log - drop_cap)];
            end
            
            % (7) Limit steepness of early decline (slope in log10 units/h)
            if any(early_win)
                dlogN_dt = gradient(log10(N + 1e-12), t);
                slope_cap = -0.20;             % don’t allow steeper than -0.20 log10/h initially
                r = [r; max(0, -(min(dlogN_dt(early_win)) - slope_cap))];
            end
        end
    end

    obj_scalar = @(x) sum(guard_resid(x).^2);

    if exist('lsqnonlin','file')==2
        opts_ls = optimoptions('lsqnonlin','Display','off', ...
            'MaxFunctionEvaluations',1500,'MaxIterations',150, ...
            'StepTolerance',1e-8,'FunctionTolerance',1e-8,'Algorithm','levenberg-marquardt','UseParallel',true);
        xhat = lsqnonlin(@(x)guard_resid(x),x0,lb,ub,opts_ls);
    elseif exist('fmincon','file')==2
        opts_fc=optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
            'MaxFunctionEvaluations',2000,'MaxIterations',200,'UseParallel',true);
        xhat=fmincon(obj_scalar,x0,[],[],[],[],lb,ub,@ecoevo_nonlincon,opts_fc);
    else
        xhat=fminsearch(obj_scalar,x0);
        xhat=min(max(xhat,lb),ub);
    end

    theta_hat = xhat;           % <-- full vector, EC50r0 estimated
    p = build_phys_from_row(row,opts);
end

function [c,ceq]=ecoevo_nonlincon(x)
    EC50s0=x(1); EC50r0=x(2); rS=x(4); rR=x(6);
    c=[EC50s0-EC50r0; rR-rS]; % ≤0 constraints
    ceq=[];
end


%% ------------------------------------------------------------------------
function [tout, yout] = simulate_pbpk_one(p, BW_kg, dose_mgkg, dose_times_h, t_end_h, y0, theta, inf_min_ami)
    if nargin<8||isempty(inf_min_ami), inf_min_ami=60; end
    p.inf_h = inf_min_ami/60;

    p.EC50s0=theta(1); p.EC50r0=theta(2); p.c=theta(3);
    p.rS=theta(4); p.Emax=theta(5); p.rR=theta(6);
    p.ku=theta(7); p.sigma=theta(8); p.d=theta(9);

    dose_times = unique(dose_times_h(dose_times_h<=t_end_h));
    pos_doses = dose_times(dose_times>=0 & dose_times<t_end_h);
    dose_mg = dose_mgkg*BW_kg;

    nD = numel(pos_doses);
    p.inf_schedule=zeros(nD,3);
    for k=1:nD
        t0=pos_doses(k);
        t1=min(t0+p.inf_h,t_end_h);
        rate=dose_mg/p.inf_h;
        p.inf_schedule(k,:)=[t0,t1,rate];
    end

    opts=odeset('RelTol',1e-6,'AbsTol',1e-9,'NonNegative',1:8);
    tBreaks=sort([pos_doses,pos_doses+p.inf_h,t_end_h]);
    tout=[];yout=[];y=y0;tPrev=0;
    for j=1:numel(tBreaks)
        tNow=tBreaks(j);
        if tNow<=tPrev,continue;end
        tspan=[tPrev,tNow];
        sol=ode45(@(t,yy)rhs_pbpk_pd_with_inf(t,yy,p),tspan,y,opts);
        npts=max(2,ceil((tspan(2)-tspan(1))*6));
        tt=linspace(tspan(1),tspan(2),npts);
        yy=deval(sol,tt).';
        tout=[tout;tt(:)];
        yout=[yout;yy];
        y=yy(end,:).';
        tPrev=tNow;
    end
    [tout,idx]=unique(tout,'stable'); yout=yout(idx,:);
end


%% ------------------------------------------------------------------------
function dy = rhs_pbpk_pd_with_inf(t,y,p)
    Cb=y(1);Amino=y(2);Ch=y(3);Ck=y(4);Cr=y(5);u=y(6);Ns=y(7);Nr=y(8);
    Rin=0;
    if isfield(p,'inf_schedule')
        for q=1:size(p.inf_schedule,1)
            if t>=p.inf_schedule(q,1) && t<p.inf_schedule(q,2)
                Rin=Rin+p.inf_schedule(q,3);
            end
        end
    end

    Sb=p.Serum_blood;He=p.Hearth;Ki=p.Kidney;Re=p.Remainder;
    Qh=p.BloodFlowRateHearth;Qk=p.BloodFlowRateKidney;Qr=p.BloodFlowRateRemainder;
    Ph=p.PartitionCoefficientHearth;Pk=p.PartitionCoefficientKidney;Pr=p.PartitionCoefficientRemainder;
    CL=p.CL;kf=p.kf;
    EC50s0=p.EC50s0;EC50r0=p.EC50r0;c=p.c;rS=p.rS;Emax=p.Emax;rR=p.rR;
    sigma=p.sigma;ku=p.ku;aRS=p.aRS;aSR=p.aSR;Kcap=p.K;d=p.d;
    Cf=p.fu*Cb;expInput=ku*u;

    dCb=(1/Sb)*(-(Qk*Cb)*Sb+(Qk*(Ck/Pk))*Ki-(Qh*Cb)+(Qh*(Ch/Ph))+ ...
        (Qr*(Cr/Pr))*Re-(Qr*Cb)*Sb+(kf*Amino)-CL*Cb);
    dAmino=Rin-p.kf*Amino;
    dCh=(1/He)*(Qh*Cb-Qh*(Ch/Ph));
    dCk=(1/Ki)*((Qk*Cb)*Sb-(Qk*(Ck/Pk))*Ki);
    dCr=(1/Re)*(-(Qr*(Cr/Pr))*Re+(Qr*Cb)*Sb);
    du=sigma*(-rR*c*exp(-c*u)*(1-(Nr+aRS*Ns)/Kcap)+ ...
        Emax*Cf*EC50r0*ku*exp(expInput)/((EC50r0*exp(expInput)+Cf)^2));
    dNs=(rS*(1-(Ns+aSR*Nr)/Kcap)-Emax*Cf/(EC50s0+Cf)-d)*Ns;
    dNr=(rR*exp(-c*u)*(1-(Nr+aRS*Ns)/Kcap)-Emax*Cf/(EC50r0*exp(expInput)+Cf)-d)*Nr;

    dy=[dCb;dAmino;dCh;dCk;dCr;du;dNs;dNr];
end


%% ------------------------------------------------------------------------
function p = build_phys_from_row(row,opts)
    if ~isfield(opts,'fu'),opts.fu=1;end
    if ~isfield(opts,'kf'),opts.kf=1;end
    if ~isfield(opts,'N_scale'),opts.N_scale=1000;end

    BW_kg=row.BW_kg;
    BV_per_kg=(row.GA<37)*0.106+(row.GA>=37)*0.080;
    p.Serum_blood=BW_kg*BV_per_kg;
    p.Hearth=0.007*BW_kg;
    p.Kidney=0.010*BW_kg;
    p.Remainder=max(BW_kg-p.Hearth-p.Kidney,0.05*BW_kg);

    CO_L_h=row.CO_L_min*60;
    frac_h=0.05;frac_k=0.20;frac_r=1-frac_h-frac_k;
    Qh=max(CO_L_h*frac_h,1e-6);
    Qk=max(CO_L_h*frac_k,1e-6);
    Qr=max(CO_L_h*frac_r,1e-6);
    sc=CO_L_h/(Qh+Qk+Qr);
    p.BloodFlowRateHearth=Qh*sc;
    p.BloodFlowRateKidney=Qk*sc;
    p.BloodFlowRateRemainder=Qr*sc;

    Vd_target=row.Vd_amikacin;
    Kp_heart=1;Kp_kidney=10;
    num=Vd_target-p.Serum_blood-Kp_heart*p.Hearth-Kp_kidney*p.Kidney;
    Kp_rem=max(num/max(p.Remainder,eps),0.05);
    p.PartitionCoefficientHearth=Kp_heart;
    p.PartitionCoefficientKidney=Kp_kidney;
    p.PartitionCoefficientRemainder=Kp_rem;

    p.CL=(row.eGFR*BW_kg)/1000*60;

    p.kf=opts.kf; p.fu=opts.fu;
    p.aRS=0.5; p.aSR=0.5;
    p.K=2.32e9*opts.N_scale;
end


function r = guard_resid(th)
    Ns0 = 4.1967e5 * opts.N_scale;
    Nr0 = 125901   * opts.N_scale;
    y0  = [0;0;0;0;0;0;Ns0;Nr0];
    N0  = Ns0 + Nr0;

    r = [];
    for g = 1:numel(opts.guard_schedules)
      G      = opts.guard_schedules{g};
      t_end  = max(opts.duration_h, 3*G.tau_h);
      times  = 0:G.tau_h:t_end;

      [t,y]  = simulate_pbpk_one(p_base,BW_kg,G.dose_mgkg,times,t_end,y0,th,G.inf_min);
      Cb=y(:,1); Amino=y(:,2); u=y(:,6); Ns=y(:,7); Nr=y(:,8); N=Ns+Nr;

      % ---------------- containment ----------------
      % (1) Peak and end below baseline (as before)
      r = [r; max(0, log10(max(N)+1e-12)   - log10(N0))];
      r = [r; max(0, log10(N(end)+1e-12)   - log10(N0))];

      % (1b) NEW: area-over-baseline (mean exceedance in log10 units)
      over     = max(0, log10(N+1e-12) - log10(N0));
      area_pen = trapz(t, over) / max(t(end)-t(1), 1e-6);
      r = [r; 3.0 * area_pen];

      % ---------------- early S survival ----------------
      early    = t <= min(12, G.tau_h);
      if any(early)
        Ns_floor = 3e6;                                    % keep S alive
        minNs    = max(min(Ns(early)), 1e-12);
        r = [r; max(0, log10(Ns_floor/minNs))];
      end

      % ---------------- on-vs-off windows ----------------
      dA  = gradient(Amino,t);
      on  = (Amino > 0 & dA >= -1e-8);                     % infusion/plateau
      off = (Amino <= 1e-6);                               % practically zero drug

      % per-capita growth (approx) from time series
      dNs = gradient(Ns,t); dNr = gradient(Nr,t);
      gS  = dNs ./ max(Ns,1e-12);                          % 1/h
      gR  = dNr ./ max(Nr,1e-12);

      % ---------------- “R killed harder early, advantage wanes with u” ----
      if any(on)
        on_idx  = find(on);
        head    = on_idx(1:min(numel(on_idx), round(6)));  % first ~6 on-samples
        tail    = on_idx(max(1, numel(on_idx)-round(6)):end);

        adv_head = mean(gR(head) - gS(head));              % negative if R killed harder
        adv_tail = mean(gR(tail) - gS(tail));              % should be less negative

        r = [r; max(0, adv_head + 0.10)];                  % demand ≥0.10 1/h more kill for R
        r = [r; max(0, -(adv_tail - adv_head))];           % advantage must WANE (increase)

        % u should rise while on-drug
        du = gradient(u,t);
        r  = [r; max(0, -mean(du(on)))];                   % mean du/dt ≥ 0 on-drug
      end

      % ---------------- off-drug dominance of S ----------------
      if any(off)
        r = [r; max(0, mean(gR(off) - gS(off)) + 0.02)];   % S grows ≥0.02 1/h faster than R off-drug
        % u should relax (not rise) off-drug
        du = gradient(u,t);
        r  = [r; max(0,  mean(du(off)))];
      end

      % ---------------- limit overly abrupt early total-kill ---------------
      if any(early)
        N_start  = max(N(find(early,1,'first')), 1e-12);
        N_min    = max(min(N(early)), 1e-12);
        drop_log = log10(N_start) - log10(N_min);
        r = [r; max(0, drop_log - 0.40)];                  % ≤0.4 log10 drop in first 12 h

        dlogN_dt = gradient(log10(N+1e-12), t);
        r = [r; max(0, -(min(dlogN_dt(early)) + 0.15))];   % slope ≥ −0.15 log10/h
      end

      % ---------------- resistant peak control (as you had) ---------------
      r = [r; max(0, log10(max(Nr)+1e-12) - log10(10*Nr0+1e-12))];
    end

    % soft biology: EC50s0 ≤ EC50r0 and rR ≤ rS (lsqnonlin has no nonlcon)
    r = [r; max(0, th(1) - th(2))];
    r = [r; max(0, th(6) - th(4))];
end
