function fit_muscle

mu = 600;           % mN/mm^3
lc0 = 2.6;          % mm
ls0 = 0.234;        % mm
lambda2 = -2.23;    % 1/mm^2
C = 2;
S = 6;
P0 = 60.86;         % mN / mm^2
k1 = 9.6;           % 1/s
k2 = 5.9;           % 1/s
k3 = 65;            % 1/s
k4 = 45;            % 1/s
M = 0.5;         % kg / mm^2
B = 0.5;           % kg / mm^2 / s
xm = 0.4;
xp = 1.33;
xmax = 1.8;
s = 0.1;
L0 = 2.7;
A = 0.125;
phi = 0.1;
T = 1;
actdur = 0.36;

datafile = 'MuscleData/sin_data.mat';
musc = load(datafile);

act = @(t) mod(t,1) < 0.36;
L = @(t) L0 + A * cos(2*pi/T * (t - phi));
V = @(t) -2*pi/T * A * sin(2*pi/T * (t - phi));

dt = 0.005;
phitest = 0:0.05:0.95;
showphi = [1 6 11 17];

pertmag = 0.1;

t = (0:0.01:3)';

if (~getvar('timedata') || ~inputyn('Use existing data?', 'default',true))
    timedata = struct([]);
    for i = 1:length(phitest)
        phi = phitest(i);
        fprintf('%d/%d.  Phi = %g\n', i,length(phitest), phi);
        
        A = 0.125/2;
        L = @(t) L0 + A * cos(2*pi/T * (t - phi));
        V = @(t) -2*pi/T * A * sin(2*pi/T * (t - phi));
        X0 = [ls0   0   0   0];
        
        %options for ode
        odeopt = odeset('RelTol',1e-6);%, 'OutputFcn',@odeplot);
        sol = ode45(@odefcn, t([1 end]), X0, odeopt);

        x1 = deval(sol,t);
        data1.t = t;
        data1.x = x1';
        data1.lc = L(t) - data1.x(:,1);
        data1.vc = V(t) - data1.x(:,2);
        data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4));
        data1.L = L;
        data1.V = V;
        timedata = makestructarray(timedata,data1);
    end
    putvar timedata;
end

figureseries('Force vs phase');
Pcall = cat(2,timedata.Pc);
plot(t, Pcall);

if (~getvar('fitdata') || ~inputyn('Use existing fit data?', 'default',true))
    phitest = musc.phi;
    t = musc.t(:);
    t = t(isfinite(t));
    fitdata = struct([]);
    A = 0.125/2;
    L = @(t) L0 - A * cos(2*pi/T * t);
    V = @(t) 2*pi/T * A * sin(2*pi/T * t);
    for i = 1:length(phitest)
        phi = phitest(i);
        fprintf('%d/%d.  Phi = %g\n', i,length(phitest), phi);

        act = @(t) ((mod(t,1) >= 0.5-phi) & (mod(t,1) < 0.5-phi+actdur));
        X0 = [ls0   0   0   0];
        
        %options for ode
        odeopt = odeset('RelTol',1e-6);%, 'OutputFcn',@odeplot);
        sol = ode45(@odefcn, [0 t(end)], X0, odeopt);

        x1 = deval(sol,t);
        data1.t = t;
        data1.x = x1';
        data1.lc = L(t) - data1.x(:,1);
        data1.vc = V(t) - data1.x(:,2);
        data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4)) / P0;
        data1.L = L;
        data1.V = V;
        data1.muscF = musc.F(:,i);
        fitdata = makestructarray(fitdata,data1);
    end
    putvar fitdata;
end

figureseries('Model and data');
Pcmod = cat(2,fitdata.Pc);
Pcdata = cat(2,fitdata.muscF);
clf;
for i = 1:size(Pcmod,2)
    addplot(fitdata(1).t+i-1, Pcmod(:,i),'r--', fitdata(1).t+i-1,Pcdata(:,i),'k-');
end

    function [hx,dhx] = h(x)
        
        exs = exp(x/s);
        hx = s * log(1 + exs);
        if (nargout == 2)
            dhx = exs ./ (1 + exs);
        end
        
    end

    function [hl,dl] = lambda(lc)
        
        l0 = 1 + lambda2 * (lc - lc0).^2;
        if (nargout == 1)
            hl = h(l0);
        else
            [hl,dhl] = h(l0);
            dl = 2.*lambda2.*(lc - lc0) .* dhl;
        end
        
    end

    function [x,dx] = xi(vc)
        
        if (nargout == 1)
            x0 = 1 + xp * h(vc) - xm * h(-vc);
            x = h(x0);
        else
            [hvcp,dhvcp] = h(vc);
            [hvcm,dhvcm] = h(-vc);
            x0 = 1 + xp * hvcp - xm * hvcm;
            [x,dhx] = h(x0);
            
            dx = (xm .* dhvcm + xp .* dhvcp) .* ...
                dhx;
        end
        
    end

    function p = Pc(lc, vc, Caf)
        
        p = P0 .* lambda(lc) .* xi(vc) .* Caf;
        
    end

    function dx = odefcn(t,x)
        
        ls = x(1,:);
        vs = x(2,:);
        Ca = x(3,:);
        Caf = x(4,:);
        
        actval = act(t);
        Lval = L(t);
        Vval = V(t);
        
        lc = Lval - ls;
        vc = Vval - vs;
        
        gact = 1./(1+exp(2-(actval-0.5)/s));
        
        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);
        dls = vs;
        
        dvs = 1/M * (Pc(lc,vc,Caf) - mu * (ls - ls0) - B * vs);
        
        dx = [dls; dvs; dCa; dCaf];
        
    end

end

  
