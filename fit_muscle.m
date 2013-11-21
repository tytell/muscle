function fit_muscle

par.T = 1;              % sec
par.mu0 = 1;
par.mu1 = 23;

par.lambda2 = -20;
par.C = 2;
par.S = 6;
par.P0 = 67;         % mN / mm^2
par.L0 = 2.94;       % mm
par.xsec = 1;           % mm^2
par.L1 = 2.682/par.L0;

par.k1 = 9;  
par.k2 = 50; 
par.k30 = 40;
par.k40 = 19.5; 
par.km1 = 5;
par.km2 = 10;
par.k5 = 100;

par.b = 10;
par.mm = 0.5;

par.alpham = 0.8;
par.alphap = 2.9;
par.alphamax = 1.8;

par.s = 0.05;

phi = 0.1;
par.actdur = 0.36;
par.A = 0.125/par.L0;

datafile = 'MuscleData/s15sines.mat';
musc = load(datafile);

par.act = @(t) mod(t,1) < par.actdur;
par.L = @(t) 1 + A1 * cos(2*pi * (t - phi));
par.V = @(t) -2*pi * A1 * sin(2*pi * (t - phi));

dt = 0.005;
phitest = 0:0.05:0.95;
showphi = [1 6 11 17];

pertmag = 0.1;

hax = gca;

mtest = [0.005 0.01 0.1 1 2 5];
btest = [0.1 0.5 1 5 10 50 100];

par.model = 'old';
par.phi = musc.phi;

par.tdat = musc.tdata;
par.Pdat = musc.Fdata;

t0 = (0:0.01:1)';
par.t = t0;

load fit_muscle.mat
if (~exist('dx0','var') || any(all(flatten(~isnan(dx0),2:4))) || ...
        inputyn('Do overview simulation again?'))
    N = length(mtest)*length(btest);
    
    [mvals,bvals] = ndgrid(mtest,btest);
    if (exist('dx','var') && exist('Pcmod','var') && ...
            all(size(dx) == [length(t0) size(mvals,1) size(mvals,2) 2]))
        dotest = squeeze(any(~isfinite(dx(:,:,1)),1));
    else
        dotest = true(size(mvals));
        dx = zeros(length(t0)*length(musc.phi),length(mtest),length(btest),3);
        Pcmod = zeros(length(t0)*length(musc.phi),length(mtest),length(btest),3);
    end
    
    [ii,jj] = find(dotest);
    dx0 = zeros(length(t0)*length(par.phi),length(ii),2);
    Pcmod0 = zeros(length(t0)*length(par.phi),length(jj),2);

    n = NaN(length(t0),length(par.phi));
    zdata = struct('x',n(:,:,[1 1 1 1 1]),'A',NaN,'L',n,'V',n,'lc',n,'vc',n,'Pc',n,'Pcdat',n);
    data0 = repmat(zdata,[1 length(ii) 2]);

    tic;
    par.model = 'lc';
    parfor k = 1:length(ii)
        i = ii(k);
        j = jj(k);

        [dx1,Pc11,data1] = fit_muscle_fcn([mtest(i); btest(j)], @muscle_ode_fcn, par);
        dx0(:,k,1) = dx1;
        Pcmod0(:,k,1) = Pc11;
        data0(1,k,1) = data1;
        
        fprintf('%d/%d (%d%%): m = %g, b = %g\n', k,N, round(k/N*100), mtest(i),btest(j));
        fprintf('   --> lsq = ');
        fprintf('%g ', nansum(dx0(:,k,1).^2,1));
        fprintf('\n');
    end
    toc;
    
    save fit_muscle.mat dx0 Pcmod0 data0;
    
    tic;
    par.model = 'ls';
    parfor k = 1:length(ii)
        i = ii(k);
        j = jj(k);

        [dx1,Pc11,data1] = fit_muscle_fcn([mtest(i); btest(j)], @muscle_ode_fcn, par);
        dx0(:,k,2) = dx1;
        Pcmod0(:,k,2) = Pc11;
        data0(1,k,2) = data1;
        
        fprintf('%d/%d (%d%%): m = %g, b = %g\n', k,N, round(k/N*100), mtest(i),btest(j));
        fprintf('   --> lsq = ');
        fprintf('%g ', nansum(dx0(:,k,2).^2,1));
        fprintf('\n');
    end
    toc;
    
    save fit_muscle.mat dx0 Pcmod0 data0;
end

dx = reshape(dx0,length(t0)*length(par.phi),length(mtest),length(btest),2);
Pcmod = reshape(Pcmod0,length(t0)*length(par.phi),length(mtest),length(btest),2);
data = reshape(data0,length(mtest),length(btest),2);

if (~exist('dxold','var'))
    par.model = 'old';
    [dxold,Pcold,dataold] = fit_muscle_fcn([0;0], @muscle_ode_fcn, par);
end
save fit_muscle.mat dx0 Pcmod0 data0 dx Pcmod data dxold Pcold dataold;

lsqold = nansum(dxold.^2);

lsq = nansum(dx.^2,1);
lsq = squeeze(lsq);

figureseries('Fit model');
clf;

subplot(2,1,1);
pcolor(log10(btest),log10(mtest),lsq(:,:,1)/lsqold);
xtick(log10(btest), cellfun(@num2str,num2cell(btest),'UniformOutput',false));
ytick(log10(mtest), cellfun(@num2str,num2cell(mtest),'UniformOutput',false));
xlabel('B');
ylabel('m');
title('lc model');
caxis([0 2]);
colorbar;

subplot(2,1,2);
pcolor(log10(btest),log10(mtest),lsq(:,:,2)/lsqold);
xtick(log10(btest), cellfun(@num2str,num2cell(btest),'UniformOutput',false));
ytick(log10(mtest), cellfun(@num2str,num2cell(mtest),'UniformOutput',false));
xlabel('B');
ylabel('m');
colorbar;
title('ls model');
caxis([0 2]);

[~,ind1] = min(flatten(lsq(:,:,1)));
[~,ind2] = min(flatten(lsq(:,:,2)));

[i1,j1] = ind2sub([size(lsq,1) size(lsq,2)], ind1);
[i2,j2] = ind2sub([size(lsq,1) size(lsq,2)], ind2);
        
figureseries('Time comparison');
subplot(2,1,1);
clf;
for k = 1:10
    addplot(t0+k-1,data(i1,j1,1).Pcdat(:,k), 'k-', ...
        t0+k-1,data(i1,j1,1).Pc(:,k), 'b--',...
        t0+k-1,data(i2,j2,2).Pc(:,k), 'g--',...
        'LineWidth',2);
    addplot(t0+k-1,dataold.Pc(:,k), 'r-', ...
        musc.tmod+k-1, musc.Fmod(:,k)/par.P0, 'm-');
end
return

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
        X0 = [1+A1   0   0   0   1];
        
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

  