function [fp,per, data] = get_limit_cycle(fcn, dt, maxper, ics, varargin)

opt.method = 'shooting';
opt.tol = 1e-4;
opt.initialcycles = 0;
opt.fixedperiod = false;
opt.TolX = 1e-6;
opt.TolFun = 1e-6;

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
    [t1,x1] = ode45(fcn, tinit, ics, odeopt);
    ics = x1(end,:);
    if (any(isnan(ics)))
        plot(t1,x1);
        error('Initial cycles do not converge');
    end
end

switch opt.method
    case 'shooting'
        tspan = [0 maxper];
        
        [fp, d] = lsqnonlin(@(x) returnval(x, fcn, tspan, odeopt, opt.fixedperiod), ics', [],[], minopt);
        
        if (~opt.fixedperiod)
            ftrn = fcn(0,fp);
            odeopt = odeset(odeopt, 'Events',@(t,x) istransverse(t,x, fp,ftrn));
        end
        sol = ode45(fcn, tspan, fp, odeopt);
        
        if (~opt.fixedperiod)
            per = diff(sol.xe);
        else
            per = maxper;
        end
        dt1 = per / ceil(per/dt);
        t1 = (0:dt1:per)';
        x = deval(sol, t1);
        
        fp = fp';
        data.t = t1;
        data.x = x';
        data.per = per;
        data.fcn = fcn;
        data.sol = sol;
 end


function dx = returnval(x0, fcn, tspan, odeopt, isfixper)

if ~isfixper
    ftrn = fcn(0,x0);
    odeopt = odeset(odeopt, 'Events',@(t,x) istransverse(t,x, x0,ftrn));
end
sol = ode45(fcn, tspan, x0, odeopt);

if (isfield(sol,'ye'))
    dx = sol.ye(:,end) - x0;
else
    dx = sol.y(:,end) - x0;
end


function [value,isterminal,direction] = istransverse(t,x, v0,ftrn)

v = x(1:length(v0));
value = ftrn' * (v-v0);
direction = 1;
isterminal = 1;



