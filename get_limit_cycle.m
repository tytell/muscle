function [fp,per, data] = get_limit_cycle(fcn, dt, maxper, ics, varargin)

opt.method = 'shooting';
opt.tol = 1e-4;
opt.initialcycles = 0;
opt.fixedperiod = false;
opt.TolX = 1e-6;
opt.TolFun = 1e-6;
opt.odetype = 'direct';

opt = parsevarargin(opt,varargin, 4, 'allowunknown');

%options for optimization
minfields = optimset;
minfields = fieldnames(minfields);
[minopt,~] = getfieldsonly(opt, minfields);

%options for ode
odefields = odeset;
odefields = fieldnames(odefields);
[odeopt,~] = getfieldsonly(opt, odefields);

if (opt.initialcycles > 0)
    tinit = [0 opt.initialcycles*maxper];
    switch opt.odetype
        case 'direct',
            [t1,x1] = ode45(fcn, tinit, ics, odeopt);
        case 'implicit',
            [t1,x1] = ode15i(fcn, tinit, ics(1,:)', ics(2,:)', odeopt);
    end
    ics = x1(end,:);
    if (any(isnan(ics)))
        plot(t1,x1);
        error('Initial cycles do not converge');
    end
end

switch opt.method
    case 'shooting'
        tspan = [0 maxper];
   
        switch opt.odetype
            case 'direct',
                [fp, d] = lsqnonlin(@(x) returnval(x, fcn, tspan, odeopt, opt.fixedperiod), ics', [],[], minopt);
                ftrn = fcn(0,fp);

                if (~opt.fixedperiod)
                    odeopt = odeset(odeopt, 'Events',@(t,x) istransverse(t,x, fp,ftrn, 0.25*diff(tspan)));
                end
                sol = ode45(fcn, tspan, fp, odeopt);
            case 'implicit',
                [~,dfdxp] = opt.Jacobian(0,ics', ones(size(ics))');
                isimp = any(dfdxp - eye(length(ics)) ~= 0);
                [fp, d] = lsqnonlin(@(x) returnvalimplicit(x, fcn, tspan, odeopt, isimp, opt.fixedperiod), ...
                    ics', [],[], minopt);
        
                fpd = zeros(size(fp));
                eqn = fcn(0,fp,fpd);
                fpd(~isimp) = -eqn(~isimp);

                [fp,fpd] = decic(fcn, 0, fp,zeros(size(fp)), fpd,~isimp);        
                ftrn = fpd;
                
                if (~opt.fixedperiod)
                    odeopt = odeset(odeopt, 'Events',@(t,x,xp) istransverse(t,x, fp,ftrn, 0.25*diff(tspan)));
                end
                sol = ode15i(fcn, tspan, fp, fpd, odeopt);
        end
        
        if (~opt.fixedperiod)
            per = diff(sol.xe);
        else
            per = maxper;
        end
        dt1 = per / ceil(per/dt);
        t1 = (0:dt1:per)';
        [x,xp] = deval(sol, t1);
        
        fp = fp';
        data.t = t1;
        data.x = x';
        data.xp = xp';
        data.per = per;
        data.fcn = fcn;
        data.sol = sol;
 end


function dx = returnval(x0, fcn, tspan, odeopt, isfixper)

if ~isfixper
    ftrn = fcn(0,x0);
    odeopt = odeset(odeopt, 'Events',@(t,x) istransverse(t,x, x0,ftrn, 0.25*diff(tspan)));
end
sol = ode45(fcn, tspan, x0, odeopt);

if (isfield(sol,'ye') && ~isempty(sol.ye))
    dx = sol.ye(:,end) - x0;
else
    dx = sol.y(:,end) - x0;
end


function dx = returnvalimplicit(x0, fcn, tspan, odeopt, isimp, isfixper)

xp0 = zeros(size(x0));
eqn = fcn(0,x0,xp0);
xp0(~isimp) = -eqn(~isimp);

[x0,xp0] = decic(fcn, 0, x0,zeros(size(x0)), xp0,~isimp);

if ~isfixper
    ftrn = xp0;
    odeopt = odeset(odeopt, 'Events',@(t,x,xp) istransverse(t,x, x0,ftrn, 0.25*diff(tspan)));
end
sol = ode15i(fcn, tspan, x0, xp0, odeopt);

if isfield(sol,'ye')
    dx = sol.ye(:,end) - x0;
else
    dx = sol.y(:,end) - x0;
end

function [value,isterminal,direction] = istransverse(t,x, v0,ftrn, minper)

v = x(1:length(v0));
value = ftrn' * (v-v0);
direction = 1;
isterminal = t > minper;



