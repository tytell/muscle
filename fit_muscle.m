function fit_muscle

%optimization 2
%m = 0.056475; b = 0.283269; lc0 = 0.965423; k1 = 6.88101; k2 = 30.2333; 
%k30 = 45.6537; k40 = 17.3768; km1 = 19.414; km2 = 6.39595    ->> sum(dx^2) = 6.350151

par.T = 1;              % sec
par.mu0 = 1;
par.mu1 = 23;

par.lambda2 = -20;
par.C = 2;
par.S = 6;
par.P0 = 67;         % mN / mm^2
par.L0 = 2.94;       % mm
par.Lis = 2.7;
par.xsec = 1;           % mm^2
par.L1 = 2.7/par.L0;
par.lc0 = 0.85; %1-(par.L0-par.Lis)/par.L0; 

par.k1 = 9;  
par.k2 = 50; 
par.k30 = 40;
par.k40 = 19.4; 
par.km1 = 15;
par.km2 = 10;
par.k5 = 100;

par.b = 0.02;
par.mm = 0.005;

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

mtest = [0.005 0.01 0.02 0.05 0.1];
btest = [0 0.01 0.05 0.1 0.5 1 2];

par.phi = musc.phi;

par.tdat = musc.tdata;
par.Pdat = musc.Fdata;

t0 = (0:0.01:1)';
par.t = t0;

% test isometric first
A1 = 0.125 / par.L0;
phi1 = 0.2;
par.L = @(t) 0; % + A1 * cos(2*pi * (t - phi1));
par.V = @(t) 0; %-2*pi * A1 * sin(2*pi * (t - phi1));
%par.mu0 = par.mu0 + par.mu1;
%par.mu1 = 0;
X0 = [0   0   0   0   1];

%options for ode
par.model = 'old2';
odeopt = odeset('RelTol',1e-5); %, 'OutputFcn',@odeplot);
sol1 = ode45(@(t,x) muscle_ode_fcn(t,x,par), [0 t0(end)+1], X0, odeopt);

t2 = (0:0.01:2)';
x2 = deval(sol1,t2);

[~,lc2,vc2,Pc2] = muscle_ode_fcn(t2',x2,par);
x2 = x2';
lc2 = lc2';
vc2 = vc2';
Pc2 = Pc2';

phi1 = 0.2;
par.L = @(t) par.L1 + A1 * cos(2*pi * (t - phi1));
par.V = @(t) -2*pi * A1 * sin(2*pi * (t - phi1));

par.model = 'ls';
odeopt = odeset('RelTol',1e-5); %, 'OutputFcn',@odeplot);
solnew = ode45(@(t,x) muscle_ode_fcn(t,x,par), [0 t0(end)+1], X0, odeopt);

xnew = deval(solnew,t2);
[~,lcn,vcn,Pcn] = muscle_ode_fcn(t2',xnew,par);

par.model = 'ls2';
odeopt = odeset('RelTol',1e-5); %, 'OutputFcn',@odeplot);
solnew = ode45(@(t,x) muscle_ode_fcn(t,x,par), [0 t0(end)+1], X0, odeopt);

xnew2 = deval(solnew,t2);
[~,lcn2,vcn2,Pcn2] = muscle_ode_fcn(t2',xnew2,par);

par.L = @(t) par.L1 + A1 * cos(2*pi * (t - phi1));
par.model = 'old2';
solold = ode45(@(t,x) muscle_ode_fcn(t,x,par), [0 t0(end)+1], X0, odeopt);

xold = deval(solold,t2);
[~,lco,vco,Pco] = muscle_ode_fcn(t2',xold,par);

xnew = xnew';
lcn = lcn';
vcn = vcn';
Pcn = Pcn';

xnew2 = xnew2';
lcn2 = lcn2';
vcn2 = vcn2';
Pcn2 = Pcn2';

xold = xold';
lco = lco';
vco = vco';
Pco = Pco';

figureseries('0.2 comparison');
subplot(2,1,1);
plot(t2,Pcn2, t2,Pcn, t2,Pco, musc.tmod,musc.Fmod(:,3)/par.P0, musc.tdata,musc.Fdata(:,3)/par.P0);
legend('reallynew','new','old','Thelma','data');

subplot(2,1,2);
plot(t2,vcn2, t2,vcn, t2,vco);

load fit_muscle.mat
if (~exist('dx0','var') || any(all(flatten(~isnan(dx0),2:4))) || ...
        inputyn('Do overview simulation again?','default',false))
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

        param = [mtest(i); btest(j); par.lc0; par.k1; par.k2; par.k30; par.k40; par.km1; par.km2];

        [dx1,Pc11,data1] = fit_muscle_fcn(param, @muscle_ode_fcn, par);
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
    par.model = 'ls2';
    parfor k = 1:length(ii)
        i = ii(k);
        j = jj(k);

        param = [mtest(i); btest(j); par.lc0; par.k1; par.k2; par.k30; par.k40; par.km1; par.km2];
        [dx1,Pc11,data1] = fit_muscle_fcn(param, @muscle_ode_fcn, par);
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
    param = [0;0; par.lc0; par.k1; par.k2; par.k30; par.k40; par.km1; par.km2];
    [dxold,Pcold,dataold] = fit_muscle_fcn(param, @muscle_ode_fcn, par);
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
caxis auto;
colorbar;

subplot(2,1,2);
pcolor(log10(btest),log10(mtest),lsq(:,:,2)/lsqold);
xtick(log10(btest), cellfun(@num2str,num2cell(btest),'UniformOutput',false));
ytick(log10(mtest), cellfun(@num2str,num2cell(mtest),'UniformOutput',false));
xlabel('B');
ylabel('m');
colorbar;
title('ls2 model');
caxis auto;

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

%one good optimum:
%m = 0.061635; b = 0.211072; lc0 = 0.969523    ->> sum(dx^2) = 0.001688
%(optimizing just the peak values)

param0 = [0.061635; 0.211072; 0.969523; par.k1; par.k2;...
          par.k30; par.k40; par.km1; par.km2];
if inputyn('Run optimization again?', 'default',false)
    i2 = 2;
    j2 = 4;
    fprintf('Start at m = %f, b = %f\n', param0(1),param0(2));
    optopt = optimset('Display','iter-detailed','FunValCheck','on', ...
        'UseParallel','always', 'DiffMinChange',0.1);
    %par.phi = par.phi(3);
    par.model = 'ls2';

    [param,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@(p) fit_muscle_fcn(p,@muscle_ode_fcn,par), param0, [], [], optopt);

    disp(param);
    save('fit_muscle_opt.mat','param','resnorm','residual','exitflag','output');
end

param(1,:) = param0;
param(1,1:3) = [0.061635 0.211072 0.969523];
param(2,:) = [0.0542    0.2802    0.9678    6.7281   23.2794   51.3537  ...
    19.3801   17.5804 6.0156];
param(3,:) = [0.0565    0.2833    0.9654    6.8810   30.2333   45.6537  ...
    17.3768   19.4140 6.3959];

par.phi = musc.phi;
[dx1,~,datafit1] = fit_muscle_fcn(param(1,:), @muscle_ode_fcn, par);
[dx2,~,datafit2] = fit_muscle_fcn(param(2,:), @muscle_ode_fcn, par);
[dx3,~,datafit3] = fit_muscle_fcn(param(3,:), @muscle_ode_fcn, par);

dx1 = reshape(dx1,size(datafit1.Pc));
dx2 = reshape(dx2,size(datafit1.Pc));
dx3 = reshape(dx3,size(datafit1.Pc));

datafit0.Pc = zeros(size(datafit1.Pc));
datafit0.Pcdat = zeros(size(datafit1.Pc));
dx0 = zeros(size(dx1));
for i = 1:length(par.phi)
    good = isfinite(musc.tmod);
    Pc1 = interp1(musc.tmod(good),musc.Fmod(good,i), t0);
    Pcdat1 = interp1(musc.tdata,musc.Fdata(:,i), t0);
    
    datafit0.Pc(:,i) = Pc1 / par.P0;
    datafit0.Pcdat(:,i) = Pcdat1 / par.P0;
    
    dx0(:,i) = (Pc1 - Pcdat1)/par.P0;
end

dxold = reshape(dxold, size(dataold.Pc));

figureseries('Fit comparison');
clf;
h(1) = subplot(2,1,1);
h(2) = subplot(2,1,2);
for k = 1:10
    addplot(h(1), t0+k-1,datafit1.Pcdat(:,k), 'k-', ...
        t0+k-1,datafit1.Pc(:,k), 'b--',...
        t0+k-1,datafit2.Pc(:,k), 'g--',...
        t0+k-1,datafit3.Pc(:,k), 'r--',...
        'LineWidth',2);
    addplot(h(1), t0+k-1,datafit0.Pc(:,k), 'm-', ...
        t0+k-1,dataold.Pc(:,k));
    
    addplot(h(2),t0+k-1,dx1(:,k), 'b--',...
        t0+k-1,dx2(:,k), 'g--',...
        t0+k-1,dx3(:,k), 'r--',...
        'LineWidth',2);
    addplot(h(2), t0+k-1,dx0(:,k), 'm-', ...
        t0+k-1,dxold(:,k), 'c-');
end

legend(h(1), 'data','1','2','3','T0','T1');
xlabel(h(2), 'Time (sec)');
ylabel(h(1), 'Pc (scaled)');
ylabel(h(2), 'Difference from data');

fprintf('Sums of squares:\n');
fprintf('My models: %g %g %g\n', sum(dx1(:).^2),sum(dx2(:).^2),sum(dx3(:).^2));
fprintf('Thelma''s reported data: %g\n', nansum(dx0(:).^2));
fprintf('My implementation of Thelma''s model: %g\n', nansum(dxold(:).^2));
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

  
