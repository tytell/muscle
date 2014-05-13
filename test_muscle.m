function test_muscle

%optimized parameters:
%2: m = 0.0542; b = 0.2802; lc0 = 0.9678; k1 = 6.7281; k2 = 23.2794; k30 = 51.3537; k40 = 19.3801; km1 = 17.5804; km2 = 6.0156    ->> sum(dx^2) = 6.056118

filename = '~etytel01/Matlab/muscle/test_muscle.mat';
nlfilename = '~etytel01/Matlab/muscle/test_muscle_nonlin.mat';
quiet = true;

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
par.actdur = 0.36;

par.act = @(t) mod(t,1) < par.actdur;
par.L = @(t) par.L1 + par.A * cos(2*pi/T * (t - par.phi));

dt = 0.005;
phitest = 0:0.05:0.95;
showphi = [1 6 11 17];

pertmag = 0.1;
    
if (~getvar('-file',filename,'data') || (~quiet && ~inputyn('Use existing data?', 'default',true)))
    n = par.T / dt + 1;
    z = zeros(n,1);
    z5 = zeros(n,5);
    
    data = struct('t',z, 'x',z, 'xp',z, 'sol',struct([]), ...
        'per',par.T, 'fcn',@(t,x) odefcn(t,x,par), ...
        'phi',0,'L',[],'V',[], 'lc',z, 'vc',z, 'Pc',z, ...
        'jfcn',@(t,x) jfcn(t,x,par), 'fx',zeros(n,5,5), 'fexp',zeros(5,1), ...
        'fmode',zeros(n,5,5));
    data = repmat(data,[1 length(phitest)]);
    
    progress(0,length(phitest),'**** Phi test');
    for i = 1:length(phitest)
        phi1 = phitest(i);
        
        par.L = @(t) par.L1 + par.A * cos(2*pi/par.T * (t - phi1));
        par.V = @(t) -2*pi/par.T * par.A * sin(2*pi/par.T * (t - phi1));
        X0 = [0   0   0   0   1];
        
        [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), dt, par.T, X0, ...
            'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
        data1.phi = phi1;
        data1.L = par.L;
        data1.V = par.V;
        data1.lc = data1.L(data1.t) - data1.x(:,1);
        data1.vc = data1.V(data1.t) - data1.x(:,2);
        data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4), par);
        
        data1 = get_floquet(data1,@(t,x) jfcn(t,x,par), 100);
        data(i) = data1;
        
        progress(i);
    end
    putvar('-file',filename,'data');
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
    
    fill([0 par.actdur par.actdur 0],[0 0 maxforce maxforce],[0.8 0.8 0.8], 'Parent',hax(i,1), ...
        'EdgeColor','none');

    fill([0 par.actdur par.actdur 0],par.L1 + par.A*[-1 -1 1 1],[0.8 0.8 0.8], 'Parent',hax(i,2), ...
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
set(gcf,'Color','w');
print('-dpdf','PhaseEffect.pdf');

fxx = cat(4, data.fx);
sgn1 = sign(fxx(1,1,1,:));
fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

fexp = cat(2,data.fexp);

figureseries('Mode time constants');
plot(phitest, log(0.5) ./ real(fexp)');
xlabel('Phase');
ylabel('t_{1/2} (sec)');
title('Mode time constants');
print('-dpdf','ModeTimeConstants.pdf');

t = data(1).t;

figureseries('Floquet modes vs. phi');
clf;
for i = 1:4,
    subplot(2,2,i);
    j = showphi(i);
    plot(data(j).t, real(fxx(:,:,1,j)));
    if any(~isreal(fxx(:,:,1,j)))
        addplot(data(j).t, imag(fxx(:,:,1,j)));
    end
    xlabel('Time (s)');
    title(sprintf('\\phi_{act} = %g',phitest(j)));
end
legend('lc','vc','Ca','Caf','m','Location','best');
print('-dpdf','FloquetModesVsPhi.pdf');

figureseries('Floquet exp');
clf;
i = 1;
subplot(2,2,1);
plot(data(i).t, fxx(:,:,1,i), 'LineWidth',2);
xlabel('Time (s)');
labellines({'\delta l_c','\delta v_c', '\delta Ca', '\delta Caf', '\delta m'}, ...
    'location',[0.05 0.1 0.6 0.5 0.8],'rotation',0);

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
print('-dpdf','FloquetExp.pdf');



figureseries('Deviation from steady vs. phi');
clf;
for i = 1:4,
    subplot(2,2,i);
    j = showphi(i);
    h1 = plot(data(j).t, data(j).x, 'k--');
    h2 = addplot(data(j).t, data(j).x + 0.1*real(fxx(:,:,1,j)));
    xlabel('Time (s)');
    title(sprintf('\\phi_{act} = %g',phitest(j)));
end
legend([h1(1); h2], 'steady','lc','vc','Ca','Caf','m','Location','best');
print('-dpdf','DevVsPhi.pdf');

if (~getvar('-file',filename,'pertdata') || (~quiet &&~inputyn('Use existing data?', 'default',true)))
    i = 3;
    t0 = data(i).t;
    xbase = data(i).x;
    Pcbase = data(i).Pc;
    lcbase = data(i).lc;
    vcbase = data(i).vc;
    
    phi1 = phitest(i);

    par.L = data(i).L;
    par.V = data(i).V;
    
    phipert = [0.2 0.7 0.2 0.2];
    pertval = [0.1 0 0 0 0; ...
        0.1 0 0 0 0;
        0 0.5 0 0 0;
        0 0 0 -0.1 0];

    progress(0,length(phipert),'**** Perturbation test');
    pertdata = struct([]);
    odeopt = odeset('RelTol',1e-6);
    for j = 1:length(phipert)
        a = find(t >= phipert(j),1);

        xinit = xbase(a,:) + pertval(j,:);
        sol = ode45(@(t,x) odefcn(t,x,par), t0(a) + [0 par.T], xinit', odeopt);

        t1 = t0;
        t1(1:a-1) = t0(1:a-1)+par.T;
        
        x1 = NaN(size(xbase));
        x1(a:end,:) = deval(sol, t0(a:end))';
        if (a > 1)
            x1(1:a-1,:) = deval(sol, par.T + t0(1:a-1))';
        end
        lc1 = par.L(t1) - x1(:,1);
        vc1 = par.V(t1) - x1(:,2);
        Pc1 = Pc(lc1, vc1, x1(:,4), par);

        pertdata(i,j).phipert = phipert(j);
        pertdata(i,j).xbase = xbase;
        pertdata(i,j).Pcbase = Pcbase;
        pertdata(i,j).lcbase = lcbase;
        pertdata(i,j).vcbase = vcbase;
        pertdata(i,j).t = t0;
        pertdata(i,j).x = x1;
        pertdata(i,j).xinit = xinit;
        pertdata(i,j).lc = lc1;
        pertdata(i,j).vc = vc1;
        pertdata(i,j).Pc = Pc1;
        progress(j);
    end
    
    putvar('-file',filename,'pertdata');
end

figureseries('Random perturbations');
clf;
i = 3;
t0 = pertdata(i,1).t;
xbase = pertdata(i,1).xbase;
Pcbase = pertdata(i,1).Pcbase;
lcbase = pertdata(i,1).lcbase;
vcbase = pertdata(i,1).vcbase;

hax = zeros(4,1);
for j = 1:4
    hax(j) = subplot(4,1,j);
end

plot(hax(1), t0,lcbase, 'k-', 'LineWidth',2);
plot(hax(2), t0,vcbase, 'k-', 'LineWidth',2);
plot(hax(3), t0,xbase(:,4), 'k-', 'LineWidth',2);
plot(hax(4), t0,Pcbase, 'k-', 'LineWidth',2);

col = 'rgbc';
for j = 1:size(pertdata,2)
    addplot(hax(1), t0,pertdata(i,j).lc, [col(j) '--']);
    addplot(hax(2), t0,pertdata(i,j).vc, [col(j) '--']);
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
print('-dpdf','Pert.pdf');

if (~getvar('-file',filename,'devdata','Pcdevall','W0','Wdev') || ...
    (~quiet && ~inputyn('Use existing data?', 'default',true)))
    devdata = struct([]);
    
    figureseries('Test');
    clf;
    
    Pcdevall = zeros(length(t),length(phitest),length(phitest));
    W0 = zeros(length(phitest),1);
    Wdev = zeros(length(phitest),length(phitest));
    
    n = 1;
    N = length(data) * length(phitest);
    progress(0,N, '**** Computing deviations...');
    odeopt = odeset('RelTol',1e-6);
    for i = 1:length(data)
        t0 = data(i).t;
        xbase = data(i).x;
        Pcbase = data(i).Pc;
        lcbase = data(i).lc;
        vcbase = data(i).vc;

        par.L = data(i).L;
        par.V = data(i).V;
        
        W0(i) = trapz(-data(i).L(t0), Pcbase);
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
            Pcdevall(:,i,j) = Pc(par.L(t0)-xfx1(:,1),par.V(t0)-xfx1(:,2),xfx1(:,4), par);

            Wdev(i,j) = trapz(-data(i).L(t), Pcdevall(:,i,j));
            %addplot(-data(i).L(t), Pcdevall(:,i,j), 'r-');
            %drawnow;
            
            if (ismember(j,showphi))
                sol = ode45(@(t,x) odefcn(t,x,par), t0(a) + [0 par.T], xfx1(a,:)', odeopt);

                x1 = NaN(size(xbase));
                t1 = t0;
                t1(1:a-1) = par.T + t0(1:a-1);
                x1(a:end,:) = deval(sol, t0(a:end))';
                if (a > 1)
                    x1(1:a-1,:) = deval(sol, par.T + t0(1:a-1))';
                end
                lc1 = par.L(t1) - x1(:,1);
                vc1 = par.V(t1) - x1(:,2);
                Pc1 = Pc(lc1,vc1, x1(:,4), par);
            else
                x1 = [];
                lc1 = [];
                vc1 = [];
                Pc1 = [];
            end

            devdata(i,j).phipert = phitest(j);
            devdata(i,j).xbase = xbase;
            devdata(i,j).lcbase = lcbase;
            devdata(i,j).vcbase = vcbase;
            devdata(i,j).Pcbase = Pcbase;
            devdata(i,j).t = t0;
            devdata(i,j).x = x1;
            devdata(i,j).xfx = xfx1;
            devdata(i,j).lc = lc1;
            devdata(i,j).vc = vc1;
            devdata(i,j).Pc = Pc1;
                        
            n = n+1;
            
            progress(n);
        end
        putvar('-file',filename,'devdata','Pcdevall','W0','Wdev');
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
        plot(devdata(i,j).t,devdata(i,j).xbase,'k-');
        addplot(devdata(i,j).t,real(devdata(i,j).x), '--', 'LineWidth',2);
        addplot(devdata(i,j).t,real(devdata(i,j).xfx),':', 'LineWidth',2);
        k = k+1;
    end
end
print('-dpdf','CheckFloquet.pdf');

figureseries('Effect of perturbations');
clf;
showphi = 1;
plot(data(showphi).t, data(showphi).Pc, 'k--','LineWidth',2);
for i = 1:4:length(phitest),
    addplot(data(showphi).t, Pcdevall(:,showphi,i), 'r-');
end
print('-dpdf','EffectPert.pdf');

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
print('-dpdf','ChangeInWork.pdf');

figureseries('Change in work contour');
clf;
contourf(phitest, phitest, bsxfun(@minus, Wdev, W0)' / max(abs(W0)));
hcol = colorbar;

xlabel('Phase of activation');
ylabel('Phase of perturbation');

ylabel(hcol, 'Fractional change');
print('-dpdf','ChangeInWorkCont.pdf');


L1test = par.lc0 + [-2*par.A 0 2*par.A];
phitest2 = 0:0.2:0.8;
if (~getvar('-file',filename,'L1data') || (~quiet && ~inputyn('Use existing data?', 'default',true)))
    L1orig = par.L1;
    L1data = struct([]);
    N = length(phitest2)*length(L1test);
    n = 1;
    progress(0,N, '**** Length test');
    for i = 1:length(phitest2)
        phi = phitest2(i);
        for j = 1:length(L1test)
            par.L1 = L1test(j);

            par.L = @(t) par.L1 + par.A * cos(2*pi/par.T * (t - phi));
            par.V = @(t) -2*pi/par.T * par.A * sin(2*pi/par.T * (t - phi));
            X0 = [0   0   0   0   1];

            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), dt, par.T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
            
            data1.lc = par.L(data1.t) - data1.x(:,1);
            data1.vc = par.V(data1.t) - data1.x(:,2);
            data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4), par);

            data1 = get_floquet(data1,@(t,x) jfcn(t,x,par), 100);
            data1.L = par.L;
            data1.V = par.V;
            L1data = makestructarray(L1data,data1);
            n = n+1;
            progress(n);
        end
    end
    L1data = reshape(L1data,[length(L1test) length(phitest2)]);
    putvar('-file',filename,'L1data');
    
    %reset L1
    par.L1 = L1orig;
end

lcall = cat(2,L1data.lc);
lcall = reshape(lcall,[],3,length(phitest2));

figureseries('Length effect');
clf;
subplot(2,2,1);
lc1 = 0.7*par.lc0:0.01:1.3*par.lc0;
plot(lc1, lambda(lc1, par), 'k--');
addplot(squeeze(lcall(:,1,:)), squeeze(lambda(lcall(:,1,:),par)),'b-', ...
    squeeze(lcall(:,2,:)), squeeze(lambda(lcall(:,2,:),par)),'g-', ...
    squeeze(lcall(:,3,:)), squeeze(lambda(lcall(:,3,:),par)),'r-', ...
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
fx = cat(3,L1data.fexp);
fx = reshape(fx,[5 size(L1data)]);
plot(phitest2,squeeze(fx(1,:,:)),'o-');
xlabel('Activation phase');
ylabel('Mode 1 exponent');
print('-dpdf','LengthEffect.pdf');

vals = fullfact([2 2 2 3]);
islen = vals(:,1) == 2;
isvel = vals(:,2) == 2;
iswork = vals(:,3) == 2;
stiffval = vals(:,4);
N = length(islen)*length(phitest2);
nldone = false;
if (getvar('-file',nlfilename','NLdata') && ...
    (numel(NLdata) == N))
    nldone = true;
end

if (~nldone || (~quiet && ~inputyn('Use existing data?', 'default',true)))
    par0 = par;
    
    progress(0,N, '**** Nonlinear calculations');
    if ~exist('NLdata','var')
        NLdata = struct([]);
        n = 0;
        i0 = 1;
    else
        i0 = floor(length(NLdata)/length(phitest2));
        n = i0*length(phitest);
        NLdata = NLdata(1:n);
        i0 = i0+1;
    end
    for i = i0:length(islen)
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
        
        for k = 1:length(phitest2)
            phi = phitest2(k);

            par.L = @(t) par.L1 + par.A * cos(2*pi/par.T * (t - phi));
            par.V = @(t) -2*pi/par.T * par.A * sin(2*pi/par.T * (t - phi));
            X0 = [0   0   0   0   1];

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
        putvar('-file',nlfilename,'NLdata');
    end
    NLdata = reshape(NLdata,[length(phitest2) length(islen)]);

    putvar('-file',filename,'NLdata');
    par = par0;
end

NLdata = reshape(NLdata,[length(phitest2) length(islen)]);

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
print('-dpdf','FloquetExpVsNonLin.pdf');


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
print('-dpdf','NonLin.pdf');

Btest = [0.05 0.1 0.28 0.5 1 2];
if (~getvar('-file',filename,'Bdata') || (~quiet && ~inputyn('Use existing B data?', 'default',true)))
    Bdata = struct([]);
    for i = 1:length(phitest2)
        phi = phitest2(i);
        for j = 1:length(Btest)
            par.b = Btest(j);
            
            par.L = @(t) par.L1 + par.A * cos(2*pi/par.T * (t - phi));
            par.V = @(t) -2*pi/par.T * par.A * sin(2*pi/par.T * (t - phi));
            X0 = [0   0   0   0   1];

            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x, par), 0.005, par.T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);

            data1.lc = par.L(data1.t) - data1.x(:,1);
            data1.vc = par.V(data1.t) - data1.x(:,2);
            data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4), par);

            data1 = get_floquet(data1,@(t,x) jfcn(t,x, par), 100);
            
            data1.L = par.L;
            data1.b = par.b;
            Bdata = makestructarray(Bdata,data1);
        end
    end
    Bdata = reshape(Bdata,[length(Btest) length(phitest2)]);
    putvar('-file',filename,'Bdata');
end

figureseries('Damping effect');
fx = cat(3,Bdata.fexp);
fx = reshape(fx,[5 size(Bdata)]);
plot(Btest,log(0.5) ./ real(squeeze(fx(1,:,:))),'o-');

xlabel('Damping coefficient');
ylabel('t_{1/2} (sec)');
title('Mode one time constant vs damping');
print('-dpdf','Damping.pdf');

k3test = 0.6:0.1:1.4;
k4test = [0.8 1 1.2];
if (~getvar('-file',filename,'k34data') || (~quiet && ~inputyn('Use existing k3 k4 data?', 'default',true)))
    phi = 0;
    k30 = par.k30;
    k40 = par.k40;
    
    par.L = @(t) par.L1 + par.A * cos(2*pi/par.T * (t - phi));
    par.V = @(t) -2*pi/par.T * par.A * sin(2*pi/par.T * (t - phi));
    X0 = [0   0   0   0   1];
    
    k34data = struct([]);
    for i = 1:length(k3test)
        par.k30 = k30*k3test(i);
        for j = 1:length(k4test)
            par.k40 = k40*k4test(j);
            
            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x, par), 0.005, par.T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);

            data1.lc = par.L(data1.t) - data1.x(:,1);
            data1.vc = par.V(data1.t) - data1.x(:,2);
            data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,4), par);

            data1 = get_floquet(data1,@(t,x) jfcn(t,x, par), 100);
            data1.L = par.L;
            data1.k30 = par.k30;
            data1.k40 = par.k40;
            k34data = makestructarray(k34data,data1);
        end
    end
    k34data = reshape(k34data,[length(k4test) length(k3test)]);
    putvar('-file',filename,'k34data');
    
    par.k30 = k30;
    par.k40 = k40;
end

figureseries('Calcium dynamics effect');
fx = cat(3,k34data.fexp);
fx = reshape(fx,[5 size(k34data)]);
plot(k3test*par.k30,log(0.5) ./ real(squeeze(fx(1,:,:))),'o-');

lab = cell(size(k4test));
for i = 1:length(k4test)
    lab{i} = sprintf('%.2g',k4test(i)*par.k40);
    if (i == 1)
        lab{i} = ['k4 = ' lab{i}];
    end
end
labellines(lab ,'rotation',0);

xlabel('k3');
ylabel('t_{1/2} (sec)');
title('Mode one time constant vs k3 and k4');
print('-dpdf','Calcium.pdf');


%--------------------------------------------------------
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

actval = par.act(t);
Lval = par.L(t);
Vval = par.V(t);

lc = Lval - ls;
vc = Vval - vs;

gact = g(actval-0.5, par);

Pcval = Pc(lc,vc,Caf, par);

dm = par.km1*Pcval.*h(-vc, par) - par.km2*(m-1).*g(vc, par);

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


    
function J = jfcn(t,x, par)

% From Mathematica:
% [0,1,0,0,0;mm.^(-1).*((-1).*mu0+(-1).*Caf.*mu1+(-1).*Caf.*alpha(V+ ...
%   (-1).*vs).*Derivative(1)(lambda)(L+(-1).*ls)),mm.^(-1).*((-1).*b+( ...
%   -1).*Caf.*lambda(L+(-1).*ls).*Derivative(1)(alpha)(V+(-1).*vs)),0, ...
%   mm.^(-1).*((-1).*ls.*mu1+alpha(V+(-1).*vs).*lambda(L+(-1).*ls)),0; ...
%   0,0,(-1).*(1+(-1).*Caf).*k30.*m.^(-1/2)+(-1).*Ca.*k2.*(1+(-1).*g(( ...
%   -0.5E0)+a))+k2.*(C+(-1).*Ca+(-1).*Caf+(-1).*S).*(1+(-1).*g(( ...
% -0.5E0)+a))+(-1).*k1.*g((-0.5E0)+a),Ca.*k30.*m.^(-1/2)+(1+(-1).* ...
%  Caf).*k40.*m.^(1/2)+(-1).*Caf.*k40.*m.^(1/2)+(-1).*Ca.*k2.*(1+(-1) ...
%   .*g((-0.5E0)+a))+(-1).*k1.*g((-0.5E0)+a),(1+(-1).*Caf).*((1/2).* ...
%   Ca.*k30.*m.^(-3/2)+(1/2).*Caf.*k40.*m.^(-1/2));0,0,(1+(-1).*Caf).* ...
%   k30.*m.^(-1/2),(-1).*Ca.*k30.*m.^(-1/2)+(-1).*(1+(-1).*Caf).*k40.* ...
%   m.^(1/2)+Caf.*k40.*m.^(1/2),(1+(-1).*Caf).*((-1/2).*Ca.*k30.*m.^( ...
%   -3/2)+(-1/2).*Caf.*k40.*m.^(-1/2));(-1).*Caf.*km1.*alpha(V+(-1).* ...
%   vs).*h((-1).*V+vs).*Derivative(1)(lambda)(L+(-1).*ls),(-1).*Caf.* ...
%   km1.*h((-1).*V+vs).*lambda(L+(-1).*ls).*Derivative(1)(alpha)(V+( ...
%   -1).*vs)+km2.*((-1)+m).*Derivative(1)(g)(V+(-1).*vs)+Caf.*km1.* ...
%   alpha(V+(-1).*vs).*lambda(L+(-1).*ls).*Derivative(1)(h)((-1).*V+ ...
%   vs),0,km1.*alpha(V+(-1).*vs).*h((-1).*V+vs).*lambda(L+(-1).*ls),( ...
%   -1).*km2.*g(V+(-1).*vs)];

ls = x(1,:);
vs = x(2,:);
Ca = x(3,:);
Caf = x(4,:);
m = x(5,:);

sqm = sqrt(m);

Lval = par.L(t);
Vval = par.V(t);

lc = Lval - ls;
vc = Vval - vs;

[lambdaval,dl] = lambda(lc, par);
[alphaval,da] = alpha(vc, par);
actval = par.act(t);
gact = g(actval-0.5, par);
[gvc,dg] = g(vc,par);

hvcm = h(-vc,par);

muval = mu(Caf, par);

J = [0,1,0,0,0; ...
    ...
    1/par.mm .* (-muval - Caf.*par.mu1 - dl.*alphaval.*Caf), ...
    1/par.mm .* (-par.b - lambdaval.*da.*Caf), ...
    0, ...
    1/par.mm .* (-ls.*par.mu1 + alphaval.*lambdaval),...
    0; ...
    ...
    0, ...
    0, ...
    -(1 - Caf).*par.k30./sqm - Ca.*par.k2.*(1 - gact) + ...
            par.k2.*(par.C - par.S - Ca - Caf) .* (1 - gact) - par.k1 .* gact, ...
    Ca.*par.k30./sqm + (1 - 2*Caf).*par.k40.*sqm - Ca.*par.k2.*(1 - gact) - par.k1.*gact, ...
    (1 - Caf) .* (0.5.*Ca .* par.k30 .* m.^-1.5 + 0.5.*Caf.*par.k40./sqm); ...
    ...
    0, ...
    0, ...
    (1 - Caf) .* par.k30./sqm, ...
    -Ca.*par.k30./sqm - (1 - 2*Caf).*par.k40.*sqm,...
    (1 - Caf) .* (-0.5.*Ca.*par.k30.*m.^-1.5 - 0.5.*Caf.*par.k40./sqm); ...
    ...
    -Caf.*par.km1.*alphaval .* hvcm .* dl, ...
    -Caf.*par.km1.*hvcm.*lambdaval.*da + par.km2.*(m-1).*dg + ...
            Caf.*par.km1.*alphaval.*lambdaval.*hvcm,...
    0,...
    par.km1.*alphaval.*hvcm.*lambdaval,...
    -par.km2.*gvc];



