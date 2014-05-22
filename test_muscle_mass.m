function test_muscle_mass

%optimized parameters:
%2: m = 0.0542; b = 0.2802; lc0 = 0.9678; k1 = 6.7281; k2 = 23.2794; k30 = 51.3537; k40 = 19.3801; km1 = 17.5804; km2 = 6.0156    ->> sum(dx^2) = 6.056118

filename = 'test_muscle_mass.mat';

quiet = true;
doplot = true;

par.L0 = 2.94;                  % mm
par.Lis = 2.7;                  % mm

par.lc0 = 1;                    % nondimensional
par.lambda2 = -20;              % nondim
par.L1 = par.lc0 - 0.05;        % nondim

par.C = 2;
par.S = 6;

par.P0 = 67;                    % mN/mm^2

par.k1 = 6.7281;                % 1/s
par.k2 = 23.2794;               % 1/s
par.k30 = 51.3537;              % 1/s
par.k40 = 19.3801;              % 1/s
par.km1 = 17.5804;              % 1/s
par.km2 = 6.0156;               % 1/s

par.mm = 0.0542;                % arbitrary
par.b = 0.2802;                 % arbitrary

par.alpham = 0.8;               
par.alphap = 2.9;
par.alphamax = 1.8;

par.mu0 = 1;
par.mu1 = 23;

par.s = 0.05;

par.A = 0.125/par.L0;
par.phi = 0.1;
par.T = 1;
par.duty = 0.36;

par.M = 0.5;
par.zeta = 0.5;
par.omegar = 2*pi*1.5;

par.act = @(t) mod(t,par.T) < par.duty;

ls1ind = 1;
vs1ind = 2;
Ca1ind = 3;
Caf1ind = 4;
m1ind = 5;
Lind = 6;
Vind = 7;

X0 = [0   0   0   0    1    ...
      0   0];
check_jacobian(0,X0', 0.02*ones(7,1), @(t,x) odefcn(t,x,par), @(t,x) jfcn(t,x,par));

if doplot
    tinit = [0 15*par.T];
    odeopt = odeset('RelTol',1e-6); %, 'OutputFcn', @odeplot);
    [t,x] = ode45(@(t,x) odefcn(t,x,par), tinit, X0, odeopt);

    lc = x(:,Lind) + par.L1 - x(:,ls1ind);
    vc = x(:,Vind) - x(:,vs1ind);

    Pcval = Pc(lc, vc, x(:,Caf1ind), par);
    figureseries('Time series');
    clf;
    xx = [0 par.duty par.duty 0]';
    yy = [0 0 1 1]';

    xx = repmat(xx,[1 15]);
    xx = xx + repmat(0:14,[4 1]);
    yy = repmat(yy,[1 15]);

    hax(1) = subplot(2,1,1);
    hold on;
    fill(xx,yy, [0.5 0.5 1]);
    plot(t,Pcval, 'LineWidth',2);
    hold off;
    ylabel('Muscle force');
    axis tight;

    hax(2) = subplot(2,1,2);
    addmplot(t,x(:,Lind),'k-', t,x(:,Vind),'k:', 'LineWidth',2);
    axis tight;

    linkaxes(hax, 'x');
    set(hax, 'XLim',[12 15]);
    set(hax(1), 'YLim',[0 1]);

    labellines(hax(2), {'L','V'});
    print('-dpdf','test_muscle_mass-0.pdf');
end

omegarvals = 2*pi* ([0.3:0.05:1 1.2 1.5 2]);
showfreq = [1 7 15 18];

if (~getvar('-file',filename,'freqdata') || (~quiet && ~inputyn('Use existing frequency data?', 'default',true)))
    freqdata = struct([]);
    
    progress(0,length(omegarvals),'**** Omegar test');
    for i = 1:length(omegarvals)
        par.omegar = omegarvals(i);
        
        X0 = [0   0   0   0    1    ...
              0   0];

        [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), 0.005, par.T, X0, ...
                'Display','final', 'fixedperiod',true, 'initialcycles',15, 'TolX',1e-8, 'RelTol',1e-6);
        data1.lc = par.L1 + data1.x(:,Lind) - data1.x(:,ls1ind);
        data1.vc = data1.x(:,Vind) - data1.x(:,vs1ind);

        data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,Caf1ind), par);

        data1.omegar = par.omegar;
        data1.zeta = par.zeta;
        data1 = get_floquet(data1,@(t,x) jfcn(t,x, par), 150);
        freqdata = makestructarray(freqdata,data1);
        
        progress(i);
    end
    
    putvar('-file',filename,'freqdata');
end

if doplot
    figureseries('Res freq');
    xx = cat(3, freqdata.x);
    Pcall = cat(2, freqdata.Pc);

    subplot(1,2,1);
    imagesc(freqdata(1).t, 1-omegarvals/(2*pi), squeeze(xx(:,Lind,:))');
    ylabel('Frequency difference f_{act} - f_{res} (Hz)');
    xlabel('Time (sec)');
    hcol = colorbar;
    ylabel(hcol, 'L');

    subplot(1,2,2);
    imagesc(freqdata(1).t, 1-omegarvals/(2*pi), Pcall');
    ylabel('Frequency difference f_{act} - f_{res} (Hz)');
    xlabel('Time (sec)');
    hcol = colorbar;
    ylabel(hcol, 'Force ratio P_c / P_0');
    print('-dpdf','test_muscle_mass-1.pdf');

    fxx = cat(4, freqdata.fx);
    sgn1 = sign(fxx(1,1,1,:));
    fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

    fexp = cat(2,freqdata.fexp);

    figureseries('Mode 1 vs res freq');
    clf;
    plot(1-omegarvals/(2*pi), log(0.5) ./ real(fexp(1,:))');
    xlabel('Frequency difference f_{act} - f_{res} (Hz)');
    ylabel('t_{1/2} (sec)');
    title('Mode one time constants');
    print('-dpdf','test_muscle_mass-2.pdf');

    figureseries('Modes vs res freq');
    clf;
    for i = 1:4,
        subplot(2,2,i);
        j = showfreq(i);

        if (isreal(fxx(:,:,1,j)))
            plot(freqdata(j).t, fxx(:,:,1,j));
        else
            plot(freqdata(j).t, real(fxx(:,:,1,j)));
            addplot(freqdata(j).t, imag(fxx(:,:,1,j)),'--');
        end
        xlabel('Time (s)');
        title(sprintf('f_{res} = %g',omegarvals(j)/(2*pi)));
    end
    legend('lc','vc','Ca','Caf','L','V','Location','best');
    print('-dpdf','test_muscle_mass-3.pdf');
end

zetaold = par.zeta;
omegarold = omegarvals;

zetavals = [0.2 1 2 4];
omegarvals = 2*pi* ([0.3 0.5 0.8 1 1.2 1.5 2]);

if (~getvar('-file',filename,'dampdata') || (~quiet && ~inputyn('Use existing damping data?', 'default',true)))
    dampdata = struct([]);
    
    X0 = [0   0   0   0    1    ...
          0   0];
    
    a = 1;
    n = length(omegarvals) * length(zetavals);
    progress(0,n,'**** Damping tests');
    for j = 1:length(zetavals)
        par.zeta = zetavals(j);
        for i = 1:length(omegarvals)
            par.omegar = omegarvals(i);
            
            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), 0.005, par.T, X0, ...
                'Display','final', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
            data1.lc = par.L1 + data1.x(:,Lind) - data1.x(:,ls1ind);
            data1.vc = data1.x(:,Vind) - data1.x(:,vs1ind);
            
            data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,Caf1ind), par);
            
            data1 = get_floquet(data1,@(t,x) jfcn(t,x,par), 150);
            data1.zeta = par.zeta;
            data1.omegar = par.omegar;
            dampdata = makestructarray(dampdata,data1);
            a = a+1;
            
            progress(a);
        end
    end
    putvar('-file',filename,'dampdata');
end

if doplot
    isomega = ismember(omegarold,omegarvals);
    islowzeta = zetavals < zetaold;

    dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)]);
    dampdata = makestructarray(dampdata(:,islowzeta), freqdata(isomega), dampdata(:,~islowzeta));
    dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)+1]);
    fexp = cat(2,dampdata.fexp);
    fexp = reshape(fexp,[size(fexp,1) length(omegarvals) length(zetavals)+1]);

    figureseries('Mode 1 vs damping');
    clf;
    plot(1-omegarvals/(2*pi), log(0.5) ./ real(squeeze(fexp(1,:,:))));
    xlabel('Frequency different f_{act} - f_{res} (Hz)');
    ylabel('t_{1/2} (sec)');
    title('Mode one time constants');

    zetavalall = [zetavals(islowzeta) zetaold zetavals(~islowzeta)];
    lab = cell(size(zetavalall));
    for i = 1:length(zetavalall),
        lab{i} = sprintf('%g',zetavalall(i));
    end
    labellines(lab,'location',[0.5 0.5 -0.2 -0.2 -0.2]);
    print('-dpdf','test_muscle_mass-4.pdf');    
end

zetavals = [0.2 1 2];
omegarvals = 2*pi* ([0.5 0.8 1 1.2 1.5 2]);
dutyvals = [0.1 0.36 0.5 0.6];
if (~getvar('-file',filename,'dutydata') || (~quiet && ~inputyn('Use existing duty cycle data?', 'default',true)))
    dutydata = struct([]);
    
    X0 = [0   0   0   0    1    ...
          0   0];
    
    a = 1;
    n = length(omegarvals) * length(zetavals) * length(dutyvals);
    progress(0,n,'**** Duty cycle tests');
    for k = 1:length(dutyvals)
        par.duty = dutyvals(k);
        for j = 1:length(zetavals)
            par.zeta = zetavals(j);
            for i = 1:length(omegarvals)
                par.omegar = omegarvals(i);

                par.act = @(t) mod(t,par.T) < par.duty;

                [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), 0.005, par.T, X0, ...
                    'Display','final', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
                data1.lc = par.L1 + data1.x(:,Lind) - data1.x(:,ls1ind);
                data1.vc = data1.x(:,Vind) - data1.x(:,vs1ind);

                data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,Caf1ind), par);

                data1 = get_floquet(data1,@(t,x) jfcn(t,x,par), 150);
                data1.zeta = par.zeta;
                data1.omegar = par.omegar;
                data1.duty = par.duty;
                dutydata = makestructarray(dutydata,data1);
                a = a+1;

                progress(a);
            end
        end
    end
    putvar('-file',filename,'dutydata');
end

vals = fullfact([2 2 2 3]);
islen = vals(:,1) == 2;
isvel = vals(:,2) == 2;
iswork = vals(:,3) == 2;
stiffval = vals(:,4);
dutyvals = [0.1 0.36 0.6];
N = length(islen) * length(dutyvals);

if (~getvar('-file',filename,'NLdata') || (~quiet && ~inputyn('Use existing data?', 'default',true)))
    par0 = par;
    
    progress(0,N, '**** Nonlinear calculations');
    NLdata = struct([]);
    n = 0;
    for i = 1:length(islen)
        if (islen(i))
            par.lambda2 = par0.lambda2;
        else
            par.lambda2 = 0;
        end
        if (isvel(i))
            par.alpham = par0.alpham;
            par.alphap = par0.alphap;
        else
            par.alpham = 0;
            par.alphap = 0;
        end
        if (iswork(i))
            par.km1 = par0.km1;
        else
            par.km1 = 0;
        end
        switch stiffval(i)
          case 1
            par.mu0 = par0.mu0 + par0.mu1;
            par.mu1 = 0;
          case 2
            par.mu0 = par0.mu0;
            par.mu1 = 0;
          case 3
            par.mu0 = par0.mu0;
            par.mu1 = par0.mu1;
        end
        
        for k = 1:length(dutyvals)
            X0 = [0   0   0   0    1    ...
                0   0];

            par.act = @(t) mod(t,par.T) < par.duty;
            
            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x, par), 0.005, par.T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);

            data1.lc = par.L(data1.t) - data1.x(:,1);
            data1.vc = par.V(data1.t) - data1.x(:,2);
            data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4), par);

            data1 = get_floquet(data1,@(t,x) jfcn(t,x, par), 100);
            data1.L = par.L;
            data1.islen = islen(i);
            data1.isvel = isvel(i);
            data1.iswork = iswork(i);
            data1.stiffval = stiffval(i);

            NLdata = makestructarray(NLdata,data1);
            n = n+1;
            progress(n);
        end
    end

    putvar('-file',filename,'NLdata');
    par = par0;
end

if doplot
    figureseries('Floquet exp vs nonlin');
    clf;
    subplot(2,3,6);
    fx = cat(3,NLdata.fexp);
    fx = reshape(fx,[5 size(NLdata)]);

    clf;
    yl = [Inf -Inf];
    hax = zeros(4,1);
    for i = 1:4,
        switch i
            case 1
                good = ~islen & ~isvel;
                ttl = 'fl=0, fv=0';
            case 2
                good = islen & ~isvel;
                ttl = 'fv=0';
            case 3
                good = ~islen & isvel;
                ttl = 'fl=0';
            case 4
                good = islen & isvel;
                ttl = 'all';
        end
        hax(i) = subplot(2,2,i);
        plot(phitest2, squeeze(real(fx(1,:,good))));
        xlabel('Activation phase');
        ylabel('Mode 1 exponent');
        title(ttl);
        if (i == 1)
            legend('none','stiff=const','work=0','both');
        end
        yl1 = ylim;
        if (yl1(1) < yl(1))
            yl(1) = yl1(1);
        end
        if (yl1(2) > yl(2))
            yl(2) = yl1(2);
        end
    end
    set(hax,'YLim',yl);
    print('-dpdf','test_muscle_mass-5.pdf');    


    Pcnl = cat(2,NLdata.Pc);
    Pcnl = reshape(Pcnl, [size(Pcnl,1) size(NLdata)]);

    figureseries('Nonlinearity effect');
    clf;
    j = 4;
    yl = [Inf -Inf];
    hax = zeros(4,1);
    for i = 1:4,
        switch i
            case 1
                good = ~islen & ~isvel;
                ttl = 'fl=0, fv=0';
            case 2
                good = islen & ~isvel;
                ttl = 'fv=0';
            case 3
                good = ~islen & isvel;
                ttl = 'fl=0';
            case 4
                good = islen & isvel;
                ttl = 'all';
        end
        hax(i) = subplot(2,2,i);
        plot(t, squeeze(Pcnl(:,j,good)));
        xlabel('Time (s)');
        ylabel('Force (mN)');
        title(ttl);
        if (i == 1)
            legend('none','stiff=const','work=0','both');
        end
        yl1 = ylim;
        if (yl1(1) < yl(1))
            yl(1) = yl1(1);
        end
        if (yl1(2) > yl(2))
            yl(2) = yl1(2);
        end
    end
    set(hax,'YLim',yl);
    print('-dpdf','test_muscle_mass-6.pdf');    
end

function [hx,dhx] = h(x, par)

exs = exp(x/par.s);
hx = par.s * log(1 + exs);
if (nargout == 2)
    dhx = exs ./ (1 + exs);
end



function [gx,dgx] = g(x, par)

e2xs = exp(-2*x/par.s);
gx = 1./(1+e2xs);

if (nargout == 2)
    dgx = -2*e2xs ./ (1+e2xs);
end



function [hl,dl] = lambda(lc, par)

l0 = 1 + par.lambda2 * (lc - par.lc0).^2;
if (nargout == 1)
    hl = h(l0, par);
else
    [hl,dhl] = h(l0, par);
    dl = 2.*par.lambda2.*(lc - par.lc0) .* dhl;
end



function [x,dx] = alpha(vc, par)

if (nargout == 1)
    x0 = 1 + par.alphap * h(vc,par) - par.alpham * h(-vc,par);
    x = h(x0,par);
else
    [hvcp,dhvcp] = h(vc,par);
    [hvcm,dhvcm] = h(-vc,par);
    x0 = 1 + par.alphap * hvcp - par.alpham * hvcm;
    [x,dhx] = h(x0,par);

    dx = (par.alpham .* dhvcm + par.alphap .* dhvcp) .* ...
        dhx;
end

function p = Pc(lc, vc, Caf, par)

p = lambda(lc,par) .* alpha(vc,par) .* Caf;


function x = mu(Caf, par)

x = par.mu0 + par.mu1*Caf;

    
function dx = odefcn(t,x, par)

ls = x(1,:);
vs = x(2,:);
Ca = x(3,:);
Caf = x(4,:);
m = x(5,:);
L = x(6,:);
V = x(7,:);

actval = par.act(t);

gact = g(actval-0.5, par);

lc = L + par.L1 - ls(1,:);
vc = V - vs(1,:);

Pcval = Pc(lc,vc,Caf, par);

dm = par.km1*Pcval .* h(-vc, par) - par.km2*(m-1).*g(vc,par);

k3 = par.k30 ./ sqrt(m);
k4 = par.k40 .* sqrt(m);

dCaf = (k3 .* Ca - k4 .* Caf) .* (1 - Caf);
dCa = (k4 .* Caf - k3 .* Ca) .* (1 - Caf) + ...
    gact .* par.k1 .* (par.C - Ca - Caf) + ...
    (1 - gact) .* par.k2 .* Ca .* (par.C - par.S - Ca - Caf);
dls = vs;

muval = mu(Caf, par);

dvs = 1/par.mm * (Pcval - par.b*vs - muval.*ls);

dL = V;
dV = 1/par.M * (-muval(1,:) .* ls(1,:) + par.b * vs(1,:) + ...
             - 2 * par.zeta * par.omegar * V - par.omegar^2 * L);

dx = [dls; dvs; dCa; dCaf; dm; dL; dV];


function J = jfcn(t,x, par)

ls = x(1,:);
vs = x(2,:);
Ca = x(3,:);
Caf = x(4,:);
m = x(5,:);
L = x(6,:);
V = x(7,:);

lc = L + par.L1 - ls;
vc = V - vs;

[lambdaval,dlambdaval] = lambda(lc,par);
[alphaval,dalphaval] = alpha(vc,par);
actval = par.act(t);

gact = g(actval-0.5, par);

[hvc,dhvc] = h(vc, par);
[gvc,dgvc] = g(vc, par);

J = zeros(7,7);

%block relating to the muscle's independent evolution
Jmusc = [ ...
    0, 1, 0, 0, 0;
    ...
    1/par.mm .* (-par.mu0 - Caf .* par.mu1 - Caf .* alphaval .* dlambdaval), ...
    1/par.mm .* (-par.b - Caf .* lambdaval .* dalphaval), ...
    0, ...
    1/par.mm .* (-par.mu1 .* ls + lambdaval .* alphaval), ...
    0; ...
    ...
    0, ...
    0, ...
    -(1 - Caf) .* par.k30 ./ sqrt(m) - Ca .* par.k2 .* (1 - gact) + ...
       par.k2 .*(par.C - par.S - Ca - Caf) .* (1 - gact) - par.k1 .* gact, ...
    Ca .* par.k30 ./ sqrt(m) + (1 - Caf).*par.k40.*sqrt(m) - Caf.*par.k40.*sqrt(m) - ...
       Ca.*par.k2 .* (1 - gact) - par.k1 .* gact, ...
    (1 - Caf) .* (Ca .* par.k30 ./ (2 * m.^1.5) + Caf .* par.k40 ./ (2*sqrt(m))); ...
    ...
    0, ...
    0, ...
    (1 - Caf) .* par.k30 ./ sqrt(m), ...
    -Ca .* par.k30 ./ sqrt(m) - (1 - Caf) .* par.k40 .* sqrt(m) + ...
       Caf .* par.k40 .* sqrt(m), ...
    (1 - Caf) .* (-Ca .* par.k30 ./ (2*m.^1.5) - Caf .* par.k40 ./ (2*sqrt(m))); ...
    ...
    -Caf .* par.km1 .* alphaval .* hvc .* dlambdaval, ...
    -Caf .* par.km1 .* hvc .* lambdaval .* dalphaval + ...
       par.km2 .* (m - 1) .* dgvc + Caf .* par.km1 .* alphaval .* lambdaval .* dhvc, ...
    0, ...
    par.km1 .* alphaval .* hvc .* lambdaval, ...
    -par.km2 .* gvc
    ];

J(1:5,1:5) = Jmusc(:,:,1);

%spring mass evolution
J(6:7,6:7) = [0, 1; ...
    -par.omegar^2/par.M, -2*par.omegar*par.zeta/par.M];

%coupling between muscle and mass
J(7,1:5) = [(-par.mu0-Caf.*par.mu1)/par.M, par.b/par.M, 0, -par.mu1.*ls/par.M, 0];
J(1:5,6:7) = [...
    0, 0; ...
    ...
    1/par.mm * Caf .* alphaval .* dlambdaval, ...
    1/par.mm * Caf .* lambdaval .* dalphaval; ...
    ...
    0, 0; ...
    0, 0; ...
    Caf .* par.km1 .* alphaval .* hvc .* dlambdaval, ...
    Caf .* par.km1 .* hvc .* lambdaval .* dalphaval - ...
       par.km2 .* (m - 1) .* dgvc - Caf .* par.km1 .* alphaval .* lambdaval .* dhvc];
    

  
