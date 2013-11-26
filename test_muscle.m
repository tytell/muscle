function test_muscle

%optimized parameters:
%2: m = 0.0542; b = 0.2802; lc0 = 0.9678; k1 = 6.7281; k2 = 23.2794; k30 = 51.3537; k40 = 19.3801; km1 = 17.5804; km2 = 6.0156    ->> sum(dx^2) = 6.056118

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

par.s = 0.05;

par.A = 0.125;
par.phi = 0.1;
par.T = 1;
par.actdur = 0.36;

par.act = @(t) mod(t,1) < par.actdur;
par.L = @(t) par.L1 + par.A * cos(2*pi/T * (t - par.phi));

dt = 0.005;
phitest = 0:0.05:0.95;
showphi = [1 6 11 17];

pertmag = 0.1;
    
if (~getvar('data') || ~inputyn('Use existing data?', 'default',true))
    data = struct([]);
    for i = 1:length(phitest)
        phi = phitest(i);
        
        L = @(t) par.L1 + par.A * cos(2*pi/par.T * (t - phi));
        X0 = [0   0   0   0   1];
        
        [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
            'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
        data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
        
        data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100);
        data1.L = L;
        data = makestructarray(data,data1);
    end
    putvar data;
end

Pcall = cat(2,data.Pc);

figureseries('Phase effect');
clf;
maxforce = max(Pcall(:));
hax = -1*ones(4,2);
for i = 1:4
    hax(i,1) = subplot(2,2,i);
    hax(i,2) = copyobj(hax(i,1), gcf);
    
    pos = get(hax(i,1),'Position');
    height = pos(4);
    pos1 = pos;
    pos1(2) = pos1(2) + 0.25*height;
    pos1(4) = 0.75*height;
    pos2 = pos;
    pos2(4) = 0.25*height;
    
    set(hax(i,1),'Position',pos1);
    set(hax(i,2),'Position',pos2);
    
    fill([0 0.36 0.36 0],[0 0 maxforce maxforce],[0.8 0.8 0.8], 'Parent',hax(i,1), ...
        'EdgeColor','none');

    fill([0 0.36 0.36 0],L0 + [-A -A A A],[0.8 0.8 0.8], 'Parent',hax(i,2), ...
        'EdgeColor','none');

    j = showphi(i);
    
    addplot(hax(i,1), data(j).t, data(j).Pc, 'r-', 'LineWidth',2);
    addplot(hax(i,2), data(j).t, data(j).L(data(j).t),'k-');
    axis(hax(i,2),'tight');
    xtick(hax(i,1), 'labeloff');
    
    xlabel(hax(i,2),'Time (s)');
    ylabel(hax(i,1),'Force (mN)');
    
    yl = get(hax(i,1),'YLim');
    maxforce = max(maxforce,yl(2));
end;
set(hax(:,1),'YLim',[0 maxforce]);
set(hax, 'TickDir','out');
set(hax(:,2), 'YAxisLocation','right');

fxx = cat(4, data.fx);
sgn1 = sign(fxx(1,1,1,:));
fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

fexp = cat(2,data.fexp);

figureseries('Mode time constants');
plot(phitest, log(0.5) ./ real(fexp)');
xlabel('Phase');
ylabel('t_{1/2} (sec)');
title('Mode time constants');

t = data(1).t;

figureseries('Floquet modes vs. phi');
clf;
for i = 1:4,
    subplot(2,2,i);
    j = showphi(i);
    plot(data(j).t, fxx(:,:,1,j));
    xlabel('Time (s)');
    title(sprintf('\\phi_{act} = %g',phitest(j)));
end
legend('lc','vc','Ca','Caf','Location','best');

figureseries('Floquet exp');
clf;
i = 1;
subplot(2,2,1);
plot(data(i).t, fxx(:,:,1,i), 'LineWidth',2);
xlabel('Time (s)');
labellines({'\delta l_c','\delta v_c', '\delta Ca', '\delta Caf'}, ...
    'location',[0.7 0.8 0.7 0.42],'rotation',0);

a = find(t >= phitest(6),1);
k = [a:length(t) 1:a-1]';
dec = exp(data(i).fexp(1)*data(i).t);
dec1 = zeros(size(dec));
dec1(k) = dec;

subplot(2,2,2);
plot(data(i).t, dec1,'k-', 'LineWidth',2);
xlabel('Time (s)');

subplot(2,4,6:7);
plot(data(i).t, bsxfun(@times, fxx(:,:,1,i), dec1), 'LineWidth',2);
xlabel('Time (s)');



figureseries('Deviation from steady vs. phi');
clf;
for i = 1:4,
    subplot(2,2,i);
    j = showphi(i);
    h1 = plot(data(j).t, data(j).x, 'k--');
    h2 = addplot(data(j).t, data(j).x + 0.1*fxx(:,:,1,j));
    xlabel('Time (s)');
    title(sprintf('\\phi_{act} = %g',phitest(j)));
end
legend([h1(1); h2], 'steady','lc','vc','Ca','Caf','Location','best');

if (~getvar('pertdata') || ~inputyn('Use existing data?', 'default',true))
    i = 3;
    t0 = data(i).t;
    xbase = data(i).x;
    Pcbase = data(i).Pc;
    
    L = data(i).L;
    
    phipert = [0.2 0.7 0.2 0.2];
    pertval = [0.1 0 0 0; ...
        0.1 0 0 0;
        0 0.5 0 0;
        0 0 0 -0.1];

    pertdata = struct([]);
    odeopt = odeset('RelTol',1e-6);
    for j = 1:length(phipert)
        a = find(t >= phipert(j),1);

        xinit = xbase(a,:) + pertval(j,:);
        sol = ode45(@odefcn, t0(a) + [0 T], xinit', odeopt);

        x1 = NaN(size(xbase));
        x1(a:end,:) = deval(sol, t0(a:end))';
        if (a > 1)
            x1(1:a-1,:) = deval(sol, T + t0(1:a-1))';
        end
        Pc1 = Pc(x1(:,1), x1(:,2), x1(:,4));

        pertdata(i,j).phipert = phipert(j);
        pertdata(i,j).x0 = xbase;
        pertdata(i,j).Pc0 = Pcbase;
        pertdata(i,j).t = t0;
        pertdata(i,j).x = x1;
        pertdata(i,j).xinit = xinit;
        pertdata(i,j).Pc = Pc1;
    end
    
    putvar pertdata;
end

figureseries('Random perturbations');
clf;
i = 3;
t0 = pertdata(i,1).t;
xbase = pertdata(i,1).x0;
Pcbase = pertdata(i,1).Pc0;

hax = tight_subplot(4,1, 0.02,[0.12 0.01],[0.12 0.01]);

plot(hax(1), t0,xbase(:,1), 'k-', 'LineWidth',2);
plot(hax(2), t0,xbase(:,2), 'k-', 'LineWidth',2);
plot(hax(3), t0,xbase(:,4), 'k-', 'LineWidth',2);
plot(hax(4), t0,Pcbase, 'k-', 'LineWidth',2);

col = 'rgbc';
for j = 1:size(pertdata,2)
    addplot(hax(1), t0,pertdata(i,j).x(:,1), [col(j) '--']);
    addplot(hax(2), t0,pertdata(i,j).x(:,2), [col(j) '--']);
    addplot(hax(3), t0,pertdata(i,j).x(:,4), [col(j) '--']);
    addplot(hax(4), t0,pertdata(i,j).Pc, [col(j) '--']);
end

linkaxes(hax(1:4), 'x');
axis(hax(1),'tight');
axis(hax(2),'tight');
axis(hax(3),'tight');
axis(hax(4),'tight');
xtick(hax(1),'labeloff');
xtick(hax(2),'labeloff');
xtick(hax(3),'labeloff');

ylabel(hax(1),'l_c (cm)');
ylabel(hax(2),'v_c (cm/s)');
ylabel(hax(3),'Caf');
ylabel(hax(4),'P_c (mN)');
xlabel(hax(4),'Time (sec)');

set(hax,'Box','off', 'TickDir','out');

if (~getvar('devdata','Pcdevall','W0','Wdev') || ~inputyn('Use existing data?', 'default',true))
    devdata = struct([]);
    
    figureseries('Test');
    clf;
    
    Pcdevall = zeros(length(t),length(phitest),length(phitest));
    W0 = zeros(length(phitest),1);
    Wdev = zeros(length(phitest),length(phitest));
    
    n = 1;
    N = length(data) * length(phitest);
    timedWaitBar(0, 'Computing deviations...');
    odeopt = odeset('RelTol',1e-6);
    for i = 1:length(data)
        t0 = data(i).t;
        xbase = data(i).x;
        Pcbase = data(i).Pc;

        L = data(i).L;
        
        W0(i) = trapz(-data(i).L(t), Pcbase);
        % plot(-data(i).L(t), Pcbase, 'k-');
        % fprintf('Total work at phase %g = %g\n', phitest(i), W0(i));

        dec = exp(data(i).fexp(1)*data(i).t);
        for j = 1:length(phitest),
            a = find(t >= phitest(j),1);
            k = [a:length(t) 1:a-1]';
            dec1 = zeros(size(dec));
            dec1(k) = dec;

            dev1 = fxx(:,:,1,i);
            %dev1 = dev1 / sqrt(sum(dev1(1,:).^2));
            dev1 = bsxfun(@times, dev1, dec1);

            xfx1 = xbase + 0.2*dev1;
            Pcdevall(:,i,j) = Pc(xfx1(:,1),xfx1(:,2),xfx1(:,4));

            Wdev(i,j) = trapz(-data(i).L(t), Pcdevall(:,i,j));
            %addplot(-data(i).L(t), Pcdevall(:,i,j), 'r-');
            %drawnow;
            
            if (ismember(j,showphi))
                sol = ode45(@odefcn, t0(a) + [0 T], xfx1(a,:)', odeopt);

                x1 = NaN(size(xbase));
                x1(a:end,:) = deval(sol, t0(a:end))';
                if (a > 1)
                    x1(1:a-1,:) = deval(sol, T + t0(1:a-1))';
                end
                Pc1 = Pc(x1(:,1), x1(:,2), x1(:,4));
            else
                x1 = [];
                Pc1 = [];
            end

            devdata(i,j).phipert = phi;
            devdata(i,j).x0 = xbase;
            devdata(i,j).Pc0 = Pcbase;
            devdata(i,j).t = t0;
            devdata(i,j).x = x1;
            devdata(i,j).xfx = xfx1;
            devdata(i,j).Pc = Pc1;
                        
            n = n+1;
            
            timedWaitBar(n/N, 'Computing deviations...');
        end
        putvar devdata Pcdevall W0 Wdev;
        %pause;
    end
end

figureseries('Check Floquet');
clf;
showphi = [1 6 11 17];
k = 1;
for ii = 1:4
    i = showphi(ii);
    for jj = 1:4
        j = showphi(jj);
        
        subplot(4,4,k);
        plot(devdata(i,j).t,devdata(i,j).x0,'k-');
        addplot(devdata(i,j).t,devdata(i,j).x, '--', 'LineWidth',2);
        addplot(devdata(i,j).t,devdata(i,j).xfx,':', 'LineWidth',2);
        k = k+1;
    end
end

figureseries('Effect of perturbations');
clf;
showphi = 1;
plot(data(showphi).t, data(showphi).Pc, 'k--','LineWidth',2);
for i = 1:4:length(phitest),
    addplot(data(showphi).t, Pcdevall(:,showphi,i), 'r-');
end

figureseries('Change in work');
clf;
showphi = 1:5:length(phitest);
plot(phitest, bsxfun(@minus,Wdev(showphi,:),W0(showphi)) / max(abs(W0)));
xlabel('Perturbation phase');
ylabel('Work');

lab = cell(size(showphi));
for i = 1:length(showphi)
    lab{i} = num2str(phitest(showphi(i)));
end

labellines(lab, 'location',[0.6 0.6 0.65 0.45],'rotation',0);

figureseries('Change in work contour');
clf;
contourf(phitest, phitest, bsxfun(@minus, Wdev, W0)' / max(abs(W0)));
hcol = colorbar;

xlabel('Phase of activation');
ylabel('Phase of perturbation');

ylabel(hcol, 'Fractional change');


L0test = lc0 + ls0 + [-2*A 0 2*A];
phitest2 = 0:0.2:0.8;
if (~getvar('L0data') || ~inputyn('Use existing data?', 'default',true))
    L0data = struct([]);
    for i = 1:length(phitest2)
        phi = phitest2(i);
        for j = 1:length(L0test)
            L0 = L0test(j);

            L = @(t) L0 + A * cos(2*pi/T * (t - phi));
            X0 = [L(0) - ls0   0   0   0];

            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
            data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));

            data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100);
            data1.L = L;
            L0data = makestructarray(L0data,data1);
        end
    end
    L0data = reshape(L0data,[length(L0test) length(phitest2)]);
    putvar L0data;
    
    %reset L0
    L0 = 2.7;
end

lcall = cat(3,L0data.x);
lcall = squeeze(lcall(:,1,:));
lcall = reshape(lcall, [size(lcall,1) size(L0data)]);

figureseries('Length effect');
clf;
subplot(2,2,1);
lc1 = 0.7*lc0:0.01:1.3*lc0;
plot(lc1, lambda(lc1), 'k--');
addplot(squeeze(lcall(:,1,:)), squeeze(lambda(lcall(:,1,:))),'b-', ...
    squeeze(lcall(:,2,:)), squeeze(lambda(lcall(:,2,:))),'g-', ...
    squeeze(lcall(:,3,:)), squeeze(lambda(lcall(:,3,:))),'r-', ...
    'LineWidth',2);
axis tight;
xlabel('lc');
ylabel('\lambda');

subplot(2,2,2);
plot(t, squeeze(lcall(:,1,:)),'b-', ...
    t, squeeze(lcall(:,2,:)),'g-', ...
    t, squeeze(lcall(:,3,:)),'r-');
xlabel('Time (s)');
ylabel('lc');

subplot(2,1,2);
fx = cat(3,L0data.fexp);
fx = reshape(fx,[4 size(L0data)]);
plot(phitest2,squeeze(fx(1,:,:)),'o-');
xlabel('Activation phase');
ylabel('Mode 1 exponent');

if (~getvar('NLdata') || ~inputyn('Use existing data?', 'default',true))
    NLdata = struct([]);
    for i = 1:2
        switch i
            case 1
                lambda2 = -2.23;
            case 2
                lambda2 = 0;
        end
        for j = 1:2
            switch j
                case 1
                    xm = 0.4;
                    xp = 1.33;
                case 2
                    xm = 0;
                    xp = 0;
            end

            for k = 1:length(phitest2)
                phi = phitest2(k);
                
                L = @(t) L0 + A * cos(2*pi/T * (t - phi));
                X0 = [L(0) - ls0   0   0   0];

                [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
                    'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
                data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));

                data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100);
                data1.L = L;

                NLdata = makestructarray(NLdata,data1);
            end
        end
    end
    NLdata = reshape(NLdata,[length(phitest2) 2 2]);
    putvar NLdata;


    lambda2 = -2.23;
    xm = 0.4;
    xp = 1.33;
end

Pcnl = cat(2,NLdata.Pc);
Pcnl = reshape(Pcnl, [size(Pcnl,1) size(NLdata)]);

figureseries('Nonlinearity effect');
clf;
for i = 1:5,
    subplot(2,3,i);
    plot(t, flatten(Pcnl(:,i,:,:),2:4));
    xlabel('Time (s)');
    ylabel('Force (mN)');
    
    if (i == 3)
        legend('normal','fv=0','fl=0','both', 'location','best');
    end
end


subplot(2,3,6);
fx = cat(3,NLdata.fexp);
fx = reshape(fx,[4 size(NLdata)]);
plot(phitest2,squeeze(fx(1,:,:)),'o-');
xlabel('Activation phase');
ylabel('Mode 1 exponent');

Btest = [5 10 15 20];
if (~getvar('Bdata') || ~inputyn('Use existing B data?', 'default',true))
    Bdata = struct([]);
    for i = 1:length(phitest2)
        phi = phitest2(i);
        for j = 1:length(Btest)
            B = Btest(j);
            
            L = @(t) L0 + A * cos(2*pi/T * (t - phi));
            X0 = [L(0) - ls0   0   0   0];

            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
            data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));

            data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100);
            data1.L = L;
            data1.B = B;
            Bdata = makestructarray(Bdata,data1);
        end
    end
    Bdata = reshape(Bdata,[length(Btest) length(phitest2)]);
    putvar Bdata;
    
    B = 10;
end

figureseries('Damping effect');
fx = cat(3,Bdata.fexp);
fx = reshape(fx,[4 size(Bdata)]);
plot(Btest,log(0.5) ./ real(squeeze(fx(1,:,:))),'o-');

xlabel('Damping coefficient');
ylabel('t_{1/2} (sec)');
title('Mode one time constant vs damping');

k3test = 0.6:0.1:1.4;
k4test = [0.8 1 1.2];
if (~getvar('k34data') || ~inputyn('Use existing k3 k4 data?', 'default',true))
    phi = 0;
    k30 = k3;
    k40 = k4;
    
    L = @(t) L0 + A * cos(2*pi/T * (t - phi));
    X0 = [L(0) - ls0   0   0   0];
    
    k34data = struct([]);
    for i = 1:length(k3test)
        k3 = k30*k3test(i);
        for j = 1:length(k4test)
            k4 = k40*k4test(j);
            
            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
            data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));

            data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100);
            data1.L = L;
            data1.k3 = k3;
            data1.k4 = k4;
            k34data = makestructarray(k34data,data1);
        end
    end
    k34data = reshape(k34data,[length(k4test) length(k3test)]);
    putvar k34data;
    
    k3 = k30;
    k4 = k40;
end

figureseries('Calcium dynamics effect');
fx = cat(3,k34data.fexp);
fx = reshape(fx,[4 size(k34data)]);
plot(k3test*k3,log(0.5) ./ real(squeeze(fx(1,:,:))),'o-');

lab = cell(size(k4test));
for i = 1:length(k4test)
    lab{i} = sprintf('%.2g',k4test(i)*k4);
    if (i == 1)
        lab{i} = ['k4 = ' lab{i}];
    end
end
labellines(lab ,'rotation',0);

xlabel('k3');
ylabel('t_{1/2} (sec)');
title('Mode one time constant vs k3 and k4');


%--------------------------------------------------------
function [hx,dhx] = h(x)

exs = exp(x/s);
hx = s * log(1 + exs);
if (nargout == 2)
    dhx = exs ./ (1 + exs);
end

function [hl,dl] = lambda(lc)

l0 = 1 + lambda2 * (lc - lc0).^2;
if (nargout == 1)
    hl = h(l0);
else
    [hl,dhl] = h(l0);
    dl = 2.*lambda2.*(lc - lc0) .* dhl;
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

function p = Pc(lc, vc, Caf)

p = P0 .* lambda(lc) .* xi(vc) .* Caf;


function dx = odefcn(t,x, par)

ls = x(1,:);
vs = x(2,:);
Ca = x(3,:);
Caf = x(4,:);
m = x(5,:);

actval = par.act(t);
Lval = par.L(t);
Vval = par.V(t);

lc = Lval - ls;
vc = Vval - vs;

gact = g(actval, par);

Pcval = Pc(lc,vc,Caf, par);

dm = par.km1*Pcval.*h(-vc, par) - par.km2*(m-1).*g(vc+0.5, par);

k3 = par.k30 ./ sqrt(m);
k4 = par.k40 .* sqrt(m);

dCaf = (k3 .* Ca - k4 .* Caf) .* (1 - Caf);
dCa = (k4 .* Caf - k3 .* Ca) .* (1 - Caf) + ...
    gact .* par.k1 .* (par.C - Ca - Caf) + ...
    (1 - gact) .* par.k2 .* Ca .* (par.C - par.S - Ca - Caf);
dls = vs;

muval = mu(Caf, par);

dvs = 1/par.mm * (Pcval - par.b*vs - muval.*ls);

dx = [dls; dvs; dCa; dCaf; dm];


    
function J = jfcn(t,x)

% From Mathematica:
% J = [0,1,0,0;M.^(-1).*((-1).*\[Mu]+(-1).*Caf.*P0.*xi(vc).*d(lambda)( ...
%  lc)),M.^(-1).*((-1).*B+(-1).*Caf.*P0.*lambda(lc).*d(xi)(vc)) ...
% ,0,(-1).*M.^(-1).*P0.*lambda(lc).*xi(vc);0,0,(-1).*(1+(-1).* ...
% Caf).*k3+(-1).*Ca.*k2.*(1+(-1).*g(act(t)))+k2.*((-1).*Ca+( ...
% -1).*Caf+Cb+(-1).*Sb).*(1+(-1).*g(act(t)))+(-1).*k1.*g(act( ...
% t)),Ca.*k3+(1+(-1).*Caf).*k4+(-1).*Caf.*k4+(-1).*Ca.*k2.*(1+ ...
% (-1).*g(act(t)))+(-1).*k1.*g(act(t));0,0,(1+(-1).*Caf).*k3,( ...
% -1).*Ca.*k3+(-1).*(1+(-1).*Caf).*k4+Caf.*k4];
%
% d(xi) = (xim.*d(h)((-1).*vc)+xip.*d(h)(vc)).*d(h)(1+(-1).*xim.*h(( ...
%    -1).*vc)+xip.*h(vc));
%
% d(lambda) = 2.*lambda2.*(lc+(-1).*lc0).*d(h)(1+lambda2.*(lc+(-1).*lc0) ...
%  .^2);
%
% d(h) = exp(1).^(s.^(-1).*x).*(1+exp(1).^(s.^(-1).*x)).^(-1)
%
% d(g) = exp(1).^(2+(-1).*((-0.5E0)+act).*s.^(-1)).*(1+exp(1).^(2+( ...
% -1).*((-0.5E0)+act).*s.^(-1))).^(-2).*s.^(-1);

lc = x(1,:);
vc = x(2,:);
Ca = x(3,:);
Caf = x(4,:);

[l,dl] = lambda(lc);
[x,dx] = xi(vc);
actval = act(t);
gact = 1./(1+exp(2-(actval-0.5)/s));

J = [0,1,0,0; ...
    ...
    1/M .* (-mu - Caf .* P0 .* x .* dl), ...
    1/M .* (-B - Caf .* P0 .* l .* dx), ...
    0, ...
    -1/M .* P0 .* l .* x; ...
    ...
    0, 0, ...
    -(1 - Caf) .* k3 - Ca .* k2 .* (1 - gact) + ...
    k2 .*(-Ca - Caf + C - S) .* (1 - gact) - k1 .* gact, ...
    Ca .* k3 + (1 - Caf).*k4 - Caf.*k4 - Ca.*k2 .* (1 - gact) - ...
    k1 .* gact;...
    ...
    0, 0, (1 - Caf).*k3, -Ca.*k3 - (1 - Caf).*k4 + Caf.*k4];



