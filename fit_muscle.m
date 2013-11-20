function fit_muscle

T = 1;              % sec
mu0 = 1;
mu1 = 23;
lambda2 = -20;
C = 2;
S = 6;
P0 = 67;         % mN / mm^2
L0 = 2.94;       % mm
xsec = 1;           % mm^2
L1 = 2.682/L0;

k1 = 9;  
k2 = 50; 
k30 = 40;
k40 = 19.5; 
km1 = 5;
km2 = 10;
k5 = 100;

b = 10;
mm = 0.5;

alpham = 0.8;
alphap = 2.9;
alphamax = 1.8;

s = 0.05;

phi = 0.1;
T = 1;
actdur = 0.36;
A = 0.125;
A1 = A/L0;

datafile = 'MuscleData/s15sines.mat';
musc = load(datafile);

act = @(t) mod(t,1) < 0.36;
L = @(t) 1 + A1 * cos(2*pi * (t - phi));
V = @(t) -2*pi * A1 * sin(2*pi * (t - phi));

dt = 0.005;
phitest = 0:0.05:0.95;
showphi = [1 6 11 17];

pertmag = 0.1;

t = (0:0.01:3)';
% 
% if (~getvar('timedata') || ~inputyn('Use existing data?', 'default',true))
%     timedata = struct([]);
%     for i = 1:length(phitest)
%         phi = phitest(i);
%         fprintf('%d/%d.  Phi = %g\n', i,length(phitest), phi);
%         
%         A = 0.125/2;
%         L = @(t) L0 + A * cos(2*pi/T * (t - phi));
%         V = @(t) -2*pi/T * A * sin(2*pi/T * (t - phi));
%         X0 = [ls0   0   0   0];
%         
%         %options for ode
%         odeopt = odeset('RelTol',1e-6);%, 'OutputFcn',@odeplot);
%         sol = ode45(@odefcn, t([1 end]), X0, odeopt);
% 
%         x1 = deval(sol,t);
%         data1.t = t;
%         data1.x = x1';
%         data1.lc = L(t) - data1.x(:,1);
%         data1.vc = V(t) - data1.x(:,2);
%         data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4));
%         data1.L = L;
%         data1.V = V;
%         timedata = makestructarray(timedata,data1);
%     end
%     putvar timedata;
% end
% 
% figureseries('Force vs phase');
% Pcall = cat(2,timedata.Pc);
% plot(t, Pcall);

figureseries('Fit model');
clf;
hax = gca;

mtest = [0.005 0.01 0.1 1 2 5];
btest = [0.1 0.5 1 5 10 50 100 500];

if (~getvar('dx','Pcmod','Pcdata') || inputyn('Do overview simulation again?'))
    n = length(mtest)*length(btest);
    t0 = (0:0.01:1)';
    
    [mvals,bvals] = ndgrid(mtest,btest);
    if (exist('dx','var') && exist('Pcmod','var') && ...
            all(size(dx) == [length(t0) size(mvals,1) size(mvals,2)]))
        dotest = squeeze(any(~isfinite(dx),1));
    else
        dotest = true(size(mvals));
        dx = zeros(length(t0)*length(musc.phi),length(mtest),length(btest));
        Pcmod = zeros(length(t0)*length(musc.phi),length(mtest),length(btest));
    end
    
    [ii,jj] = find(dotest);
   
    for k = 1:length(ii)
        i = ii(k);
        j = jj(k);
        fprintf('%d/%d (%d%%): m = %g, b = %g\n', k,n, round(k/n*100), mtest(i),btest(j));
            
        [dx1,Pc11,Pcdata] = fitfcn([mtest(i); btest(j)], []);
        dx(:,i,j) = dx1;
        Pcmod(:,i,j) = Pc11;
            
        fprintf('   --> lsq = %g\n', sum(dx(:,i,j).^2,1));
        putvar dx Pcmod Pcdata;
    end
end
        
lsq = sum(dx.^2,1);
lsq = squeeze(lsq);
pcolor(log10(btest),log10(mtest),lsq);

return;

        
optopt = optimset('Display','iter-detailed','FunValCheck','on', ...
    'Algorithm','levenberg-marquardt', 'UseParallel','always');
param = lsqnonlin(@(p) fitfcn(p,[]), [mm; b], [], [], optopt);
return;

if (~getvar('fitdata') || ~inputyn('Use existing fit data?', 'default',true))
    phitest = musc.phi;
    t = musc.tmod(:);
    t = t(isfinite(t));
    fitdata = struct([]);
    for i = 1:length(phitest)
        phi = phitest(i);
        fprintf('%d/%d.  Phi = %g\n', i,length(phitest), phi);

        A1 = 0.125/2 / L0;
        L = @(t) 1 + A1 * cos(2*pi * (t - phi));
        V = @(t) -2*pi * A1 * sin(2*pi * (t - phi));
        X0 = [0   0   0   0   1];
        
        %options for ode
        odeopt = odeset('RelTol',1e-6); %, 'OutputFcn',@odeplot);
        sol1 = ode45(@odefcn1, [0 t(end)+1], X0, odeopt);
        odeopt = odeset('RelTol',1e-6); %, 'OutputFcn',@odeplot);
        sol2 = ode23(@odefcn2, [0 t(end)+1], X0, odeopt);

        x1 = deval(sol1,t+1);
        data1.t = t;
        data1.x = x1';
        data1.lc = L(t) - data1.x(:,1);
        data1.vc = V(t) - data1.x(:,2);
        data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4));
        data1.L = L;
        data1.V = V;
        data1.muscF = interp1(musc.tdata,musc.Fdata(:,i), t);
        
        good = isfinite(musc.tmod);
        data1.oldF = interp1(musc.tmod(good),musc.Fmod(good,i), t);
        
        data1.xold = deval(sol2,t+1)';
        fitdata = makestructarray(fitdata,data1);
    end
    putvar fitdata;
end

figureseries('Model and data');
Pcmod = cat(2,fitdata.Pc);
Pcold = cat(2,fitdata.oldF);
Pold = cat(2,fitdata.xold);
Pold = Pold(:,1:5:end);
Pcdata = cat(2,fitdata.muscF);
clf;
hax(1) = subplot(2,1,1);
hax(2) = subplot(2,1,2);
for i = 1:size(Pcmod,2)
    addplot(hax(1),fitdata(1).t+i-1, Pcmod(:,i),'r--', ...
        fitdata(1).t+i-1, Pold(:,i), 'b:', fitdata(1).t+i-1, Pcold(:,i)/P0, 'g:', ...
        fitdata(1).t+i-1,Pcdata(:,i)/P0,'k-', 'LineWidth',2);
    addplot(hax(2),fitdata(1).t+i-1, fitdata(i).vc,'g-');
end

    function [hx,dhx] = h(x)
        
        exs = exp(x/s);
        hx = s * log(1 + exs);
        if (nargout == 2)
            dhx = exs ./ (1 + exs);
        end
        
    end

    function gx = g(x)
        
        gx = 1./(1+exp(-2*(x-0.5)/s));
    end

    function [hl,dl] = lambda(lc)
        
        l0 = 1 + lambda2 * (lc - 1).^2;
        if (nargout == 1)
            hl = h(l0);
        else
            [hl,dhl] = h(l0);
            dl = 2.*lambda2.*(lc - 1) .* dhl;
        end
        
    end

    function [x,dx] = alpha(vc)
        
        if (nargout == 1)
            x = zeros(size(vc));
            x(vc >= 0) = 1 + alphap * vc(vc >= 0);
            x(vc < 0) = 1 + alpham * vc(vc < 0);
            x(vc > alphamax) = alphamax;
            %x0 = 1 + alphap * h(vc) - alpham * h(-vc);
            %x = h(x0);
            %x = alphamax - h(alphamax - x);
        else
            [hvcp,dhvcp] = h(vc);
            [hvcm,dhvcm] = h(-vc);
            x0 = 1 + alphap * hvcp - alpham * hvcm;
            [x,dhx] = h(x0);
            
            dx = (alpham .* dhvcm + alphap .* dhvcp) .* ...
                dhx;
        end
        
    end

    function x = mu(Caf)
        
        x = mu0 + mu1*Caf;
        
    end

    function p = Pc(lc, vc, Caf)
        
        p = lambda(lc) .* alpha(vc) .* Caf;
        
    end

    function dx = odefcn1(t,x)
        
        ls = x(1,:);
        vs = x(2,:);
        Ca = x(3,:);
        Caf = x(4,:);
        m = x(5,:);
        
        actval = act(t);
        Lval = L(t);
        Vval = V(t);
        
        lc = Lval - ls;
        vc = Vval - vs;
        
        gact = actval; %g(actval);

        Pcval = Pc(lc,vc,Caf);
        
        dm = km1*Pcval*h(-vc) - km2*(m-1)*g(vc+0.5);
        
        k3 = k30 / sqrt(m);
        k4 = k40 * sqrt(m);
        
        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);
        dls = vs;
        
        muval = mu(Caf);
        
        dvs = 1/mm * (Pcval - b*vs - muval*ls);
        
        dx = [dls; dvs; dCa; dCaf; dm];
        
    end

    function dx = odefcn2(t,x)
        
        P = x(1,:);
        Ca = x(3,:);
        Caf = x(4,:);
        m = x(5,:);
        
        actval = act(t);
        Lval = L(t);
        Vval = V(t);
        
        gact = actval; %g(actval);

        k3 = k30 / sqrt(m);
        k4 = k40 * sqrt(m);
        
        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);
        
        muval = mu(Caf);

        lc = Lval - P/muval;
        vc_sign = muval.* Vval - k5*(Caf.*lambda(lc) - P) + P.*mu1./muval.*dCaf;
        if (vc_sign < 0)
            alpha1 = alpham;
        else
            alpha1 = alphap;
        end
        
        vc = vc_sign./(muval + k5.*Caf.*lambda(lc).*alpha1);

        Pcval = Pc(lc,vc,Caf);
        dP = k5*(Pcval-P);
        
        %dP = (lambda(lc) * Caf * (1 + alpha1*Vval + alpha1*mu1*P * dCaf/muval^2) - P) / ...
        %    (1/k5 + lambda(lc) * alpha1 * Caf/muval);
                
        dm = km1*Pcval*h(-vc) - km2*(m-1)*g(vc+0.5);
        
        dx = [dP; 0; dCa; dCaf; dm];
        
    end

    function [dx,Pcmb,Pcdat] = fitfcn(param, hax)

        mm = param(1);
        b = param(2);
        
        t1 = (0:0.01:1)';
        phitest1 = musc.phi;
        
        if (~isempty(hax))
            cla(hax,'reset');
            axis([0 length(phitest1) -1 1.8]);
            title(sprintf('mm = %f, b = %f',mm,b));
        end
        dx = zeros(length(t1),length(phitest1));
        Pcmb = zeros(length(t1),length(phitest1));
        Pcdat = zeros(length(t1),length(phitest1));
        for iff = 1:length(phitest1)
            phi1 = phitest1(iff);

            A1 = 0.125 / L0;
            L = @(t) L1 + A1 * cos(2*pi * (t - phi1));
            V = @(t) -2*pi * A1 * sin(2*pi * (t - phi1));
            X0 = [0   0   0   0   1];

            %options for ode
            odeopt = odeset('RelTol',1e-5); %, 'OutputFcn',@odeplot);
            sol1 = ode23(@odefcn1, [0 t(end)+1], X0, odeopt);

            x1 = deval(sol1,t1+1);
            x1 = x1';
            lc1 = L(t1) - x1(:,1);
            vc1 = V(t1) - x1(:,2);
            Pcmb(:,iff) = Pc(lc1, vc1, x1(:,4));
            Pcdat(:,iff) = interp1(musc.tdata,musc.Fdata(:,iff), t1) / P0;
            Pcdat(isnan(Pcdat(:,iff)),iff) = 0;
            
            dx(:,iff) = Pcmb(:,iff) - Pcdat(:,iff);
            if (~isempty(hax))
                addplot(hax,t1+iff-1,Pcmb(:,iff), t1+iff-1,Pcdat(:,iff), t1+iff-1,dx(:,iff));
                drawnow;
            end
        end
        dx = dx(:);
        Pcmb = Pcmb(:);
        Pcdat = Pcdat(:);
    end        
end

  
