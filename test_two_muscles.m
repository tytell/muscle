function test_two_muscles

mu = 600;
lc0 = 2.6;
ls0 = 0.234;
lambda2 = -2.23;
C = 2;
S = 6;
P0 = 60.86;
k1 = 9.6;
k2 = 5.9;
k3 = 65;
k4 = 45;
mm = 0.05;
B = 10;
xm = 0.4;
xp = 1.33;
ximax = 1.8;
s = 0.1;
dphi = 0.5;
T = 1;
duty = 0.36;

L0 = 2.7;
M = 0.5;
zeta = 0.5;
omegar = 2*pi*1.5;

inputyn('','clearsaved');

dt = 0.005;

omegarvals = 2*pi* ([0.3:0.05:1 1.2 1.5 2]);
showfreq = [1 7 15 18];

X0 = [L0 - ls0   0   0   0    L0 - ls0  0   0   0   L0   0];
check_jacobian(0,X0', 0.02*ones(10,1), @odefcn, @jfcn);

tinit = [0 15*T];
odeopt = odeset('RelTol',1e-6); %, 'OutputFcn', @odeoutput);
[t,x] = ode45(@odefcn, tinit, X0, odeopt);

Pcl = Pc(x(:,1), x(:,2), x(:,4));
Pcr = Pc(x(:,5), x(:,6), x(:,8));
figureseries('Time series');
clf;
xx = [0 duty duty 0; 0.5 duty+0.5 duty+0.5 0.5]';
yy = [0 0 30 30; 0 0 30 30]';

xx = repmat(xx,[1 15]);
xx(:,1:2:end) = xx(:,1:2:end) + repmat(0:14,[4 1]);
xx(:,2:2:end) = xx(:,2:2:end) + repmat(0:14,[4 1]);
yy = repmat(yy,[1 15]);

hax(1) = subplot(2,1,1);
hold on;
fill(xx(:,1:2:end),yy(:,1:2:end), [0.5 0.5 1]);
fill(xx(:,2:2:end),yy(:,2:2:end), [0.5 1 0.5]);
plot(t,Pcl, t,Pcr,'LineWidth',2);
hold off;
ylabel('Left muscle');
axis tight;

hax(2) = subplot(2,1,2);
addmplot(t,x(:,9)-L0,'k-', t,x(:,10),'k:', 'LineWidth',2);
axis tight;

linkaxes(hax, 'x');
set(hax, 'XLim',[12 15]);
set(hax(1), 'YLim',[0 30]);

labellines(hax(1), {'left','right'}, 'location',[13.12 13.62]);
labellines(hax(2), {'L','V'});

if (~getvar('freqdata') || ~inputyn('Use existing data?', 'default',true))
    freqdata = struct([]);
    
    timedWaitBar(0, 'Estimating limit cycles...');
    for i = 1:length(omegarvals)
        omegar = omegarvals(i);
        
        [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',15, 'TolX',1e-8, 'RelTol',1e-6);
        timedWaitBar((2*i-1)/(2*length(omegarvals)));
        
        data1.Pc(:,1) = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
        data1.Pc(:,2) = Pc(data1.x(:,5), data1.x(:,6), data1.x(:,8));
        data1.omegar = omegar;
        data1.zeta = zeta;
        data1 = get_floquet(data1,@(t,x) jfcn(t,x), 150);
        freqdata = makestructarray(freqdata,data1);
        timedWaitBar((2*i)/(2*length(omegarvals)));
    end
    timedWaitBar(1);
    
    putvar freqdata;
end

figureseries('Time series res freq');
clf;
xx = [0 duty duty 0; 0.5 duty+0.5 duty+0.5 0.5]';
yy = [0 0 30 30; 0 0 30 30]';

xx = repmat(xx,[1 2]);
xx(:,1:2:end) = xx(:,1:2:end) + repmat(0:1,[4 1]);
xx(:,2:2:end) = xx(:,2:2:end) + repmat(0:1,[4 1]);
yy = repmat(yy,[1 2]);

ex = [5 15 17];
t = freqdata(ex(1)).t;
Pc1 = cat(2,freqdata(ex).Pc);
Pc1 = Pc1(:,1:2:end);
Pc2 = freqdata(ex(2)).Pc(:,2);

L1 = zeros(size(Pc1));
for i = 1:3
    L1(:,i) = freqdata(ex(i)).x(:,9) - L0;
end

hax(1) = subplot(2,1,1);
hold on;
fill(xx(:,1:2:end),yy(:,1:2:end), [0.7 0.7 0.7], 'EdgeColor','none');
fill(xx(:,2:2:end),yy(:,2:2:end), 'w', 'EdgeColor','k');
hln = plot([t; t+1],[Pc1; Pc1],'LineWidth',2);
hln(end+1) = addplot([t; t+1],[Pc2; Pc2], 'g--');
hold off;
axis tight;

hax(2) = subplot(2,1,2);
addplot([t; t+1],[L1; L1], 'LineWidth',2);
axis tight;

linkaxes(hax, 'x');
labellines(hln,{'2','1','f_{act}/f_{res} = 0.67','right'}, 'location',[0.08 0.25 0.35, 0.8], ...
    'rotation',0)

figureseries('Effect of resonant frequency');
xx = cat(3, freqdata.x);
Pcall = cat(3, freqdata.Pc);

subplot(1,2,1);
imagesc(freqdata(1).t, 1./(omegarvals/(2*pi)), squeeze(xx(:,9,:))');
ylabel('Frequency ratio f_{act}/f_{res} (Hz)');
xlabel('Time (sec)');
hcol = colorbar;
ylabel(hcol, 'L');

subplot(1,2,2);
imagesc(freqdata(1).t, 1./(omegarvals/(2*pi)), squeeze(Pcall(:,1,:))' / P0);
ylabel('Frequency ratio f_{act}/f_{res} (Hz)');
xlabel('Time (sec)');
hcol = colorbar;
ylabel(hcol, 'Force ratio P_{c,1} / P_0');

fxx = cat(4, freqdata.fx);
sgn1 = sign(fxx(1,1,1,:));
fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

fexp = cat(2,freqdata.fexp);

figureseries('Floquet exponents vs resonant frequency');
clf;
plot(1./(omegarvals/(2*pi)), log(0.5) ./ real(fexp(1:2,:))');
xlabel('Frequency ratio f_{act}/f_{res} (Hz)');
ylabel('t_{1/2} (sec)');
title('Mode one time constants');

figureseries('Floquet modes vs resonant frequency');
clf;
hln = -1*ones(10,1);
for i = 1:4,
    subplot(2,2,i);
    j = showfreq(i);
    
    if (isreal(fxx(:,:,1,j)))
        hln(1:4) = plot(freqdata(j).t, fxx(:,1:4,1,j));
        hln(5:8) = addplot(freqdata(j).t, fxx(:,5:8,1,j),'--');
        hln(9:10) = addmplot(freqdata(j).t, fxx(:,9:10,1,j),'k-|k--','LineWidth',2);
    else
        plot(freqdata(j).t, real(fxx(:,:,1,j)));
        addplot(freqdata(j).t, imag(fxx(:,:,1,j)),'--');
    end
    xlabel('Time (s)');
    title(sprintf('f_{res} = %g',omegarvals(j)/(2*pi)));
end
legend(hln([1:4 9 10]),'lc','vc','Ca','Caf','L','V','Location','best');

odeopt = odeset('RelTol',1e-6); %, 'OutputFcn', @odeoutput);
omegar = 2*pi*0.5;
i = find(omegarvals == omegar);
ph = 0.2;
j = first(freqdata(i).t >= ph);

X0 = freqdata(i).x(j,:);
X0 = X0 + fxx(j,:,1,i)*0.3;
[t,x] = ode45(@odefcn, freqdata(i).t(j) + [0 2], X0, odeopt);

Pc1 = freqdata(i).Pc;
Pcpert(:,1) = Pc(x(:,1), x(:,2), x(:,4));
Pcpert(:,2) = Pc(x(:,5), x(:,6), x(:,8));

figureseries('Position perturbation');
hax = zeros(2,1);
hax(1) = subplot(2,1,1);
plot([freqdata(i).t; freqdata(i).t+1], [freqdata(i).x(:,9); freqdata(i).x(:,9)]-L0, 'k-', ...
    'LineWidth',1.5);
addplot(t,x(:,9)-L0, 'LineWidth',2);

ylabel('L (cm)');
xtick labeloff;

hax(2) = subplot(2,1,2);
plot([freqdata(i).t; freqdata(i).t+1], [freqdata(i).Pc; freqdata(i).Pc], 'k-', ...
    'LineWidth',1.5);
addmplot(t,Pcpert,'r-|r--', 'LineWidth',2);

ylabel('P_c (mN)');
xlabel('Time (sec)');


zetaold = zeta;
omegarold = omegarvals;

zetavals = [0.2 1 2 4];
omegarvals = 2*pi* ([0.3 0.5 0.8 1 1.2 1.5 2]);

if (~getvar('dampdata') || ~inputyn('Use existing data?', 'default',true))
    dampdata = struct([]);
    
    X0 = [L0 - ls0   0   0   0    L0 - ls0  0   0   0   L0   0];
    
    timedWaitBar(0,'Calculating damping and resonance');
    a = 1;
    n = length(omegarvals) * length(zetavals);
    for j = 1:length(zetavals)
        zeta = zetavals(j);
        for i = 1:length(omegarvals)
            omegar = omegarvals(i);
            
            fprintf('Zeta = %f, OmegaR = %f\n', zeta, omegar);
            
            [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
            data1.Pc(:,1) = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
            data1.Pc(:,2) = Pc(data1.x(:,5), data1.x(:,6), data1.x(:,8));
            
            data1 = get_floquet(data1,@(t,x) jfcn(t,x), 150);
            data1.zeta = zeta;
            data1.omegar = omegar;
            dampdata = makestructarray(dampdata,data1);
            timedWaitBar(a/n,'Calculating damping and resonance');
            a = a+1;
            putvar dampdata;
        end
    end
end

isomega = ismember(omegarold,omegarvals);
islowzeta = zetavals < zetaold;

dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)]);
dampdata = makestructarray(dampdata(:,islowzeta), freqdata(isomega), dampdata(:,~islowzeta));
dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)+1]);
fexp = catuneven(2,dampdata.fexp);
fexp = reshape(fexp,[size(fexp,1) length(omegarvals) length(zetavals)+1]);

for i = 1:2
    switch i
        case 1
            figureseries('Damping and resonance, mode 1');
        case 2
            figureseries('Damping and resonance, mode 2');
    end            
    clf;
    hax(i) = gca;
    plot(1./(omegarvals/(2*pi)), log(0.5) ./ real(squeeze(fexp(i,:,:))));
    xlabel('Frequency ratio f_{act}/f_{res} (Hz)');
    ylabel('t_{1/2} (sec)');
    title(sprintf('Mode %d time constants',i));
end

zetavalall = [zetavals(islowzeta) zetaold zetavals(~islowzeta)];
lab = cell(size(zetavalall));
for i = 1:length(zetavalall),
    lab{i} = sprintf('%g',zetavalall(i));
end
labellines(hax(1), lab,'location',-0.2);


figureseries('Damping and floquet modes');
clf;

showdamp = [1 1; length(omegarvals) 1; 1 length(zetavalall); length(omegarvals) length(zetavalall)];
for k = 1:4
    subplot(2,2,k);
    i = showdamp(k,1);
    j = showdamp(k,2);
    
    hln = zeros(6,1);
    hln(1:4) = plot(dampdata(i,j).t, dampdata(i,j).fx(:,1:4,1));
    addplot(dampdata(i,j).t, dampdata(i,j).fx(:,5:8,1),'--');
    hln(5:6) = addmplot(dampdata(i,j).t, dampdata(i,j).fx(:,9:10,1), 'k-|k--', 'LineWidth',2);
    xlabel('Time (sec)');
    
    title(sprintf('f_{act}/f_{res} = %g, \\zeta = %g, t_{1/2} = %g', ...
        1./(omegarvals(i)/(2*pi)), zetavalall(j), log(0.5) ./ real(fexp(1,i,j))));
end
labellines(hln, {'lc','vc','Ca','Caf','L','V'});

figureseries('Effect of the mode');
clf;

i = 1;
j = length(zetavalall);
plot(dampdata(i,j).t, dampdata(i,j).fx(:,1:4,1));
addplot(dampdata(i,j).t, dampdata(i,j).fx(:,5:8,1),'--');
addmplot(dampdata(i,j).t, dampdata(i,j).fx(:,9:10,1), 'k-|k--', 'LineWidth',2);
xlabel('Time (sec)');
ylabel('Deviation');

labellines(hln, {'lc1','vc1','Ca1','Caf1','lc2','vc2','Ca2','Caf2','L','V'});

zetaoldvals = zetavals;
omegaroldvals = omegarvals;
dutycycleold = 0.36;

dutycyclevals = [0.2 0.5 0.7];
zetavals = [0.2 4];
omegarvals = 2*pi* ([0.3 2]);
if (~getvar('dutydata') || ~inputyn('Use existing data?', 'default',true))
    dutydata = struct([]);
    
    X0 = [L0 - ls0   0   0   0    L0 - ls0  0   0   0   L0   0];
    
    timedWaitBar(0,'Calculating duty cycles');
    a = 1;
    n = length(omegarvals) * length(zetavals) * length(dutycyclevals);
    for k = 1:length(dutycyclevals)
        duty = dutycyclevals(k);
        
        for j = 1:length(zetavals)
            zeta = zetavals(j);
            for i = 1:length(omegarvals)
                omegar = omegarvals(i);
                fprintf('Duty = %g, Zeta = %g, OmegaR = %g\n', duty, zeta, omegar);

                [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, X0, ...
                    'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
                data1.Pc(:,1) = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
                data1.Pc(:,2) = Pc(data1.x(:,5), data1.x(:,6), data1.x(:,8));

                data1 = get_floquet(data1,@(t,x) jfcn(t,x), 150);
                
                data1.dutycycle = duty;
                data1.zeta = zeta;
                data1.omegar = omegar;
                dutydata = makestructarray(dutydata,data1);
                timedWaitBar(a/n);
                a = a+1;
                putvar dutydata;
            end
        end
    end
end

isomega = ismember(omegaroldvals,omegarvals);
iszeta = ismember(zetaoldvals,zetavals);
islowduty = dutycyclevals < dutycycleold;

dutydata = reshape(dutydata, [length(omegarvals) length(zetavals) length(dutycyclevals)]);
dutydata = makestructarray(dutydata(:,:,islowduty), dampdata(isomega,iszeta), dutydata(:,:,~islowduty));
dutydata = reshape(dutydata, [length(omegarvals) length(zetavals) length(dutycyclevals)+1]);

dutyall = [dutycyclevals(islowduty) dutycycleold dutycyclevals(~islowduty)];

figureseries('Force length vs duty cycle');
xx = cat(3, dutydata.x);
Lall = squeeze(xx(:,9,:));
Lall = reshape(Lall,[size(Lall,1) 4 4]);

Pcall = cat(3, dutydata.Pc);
Pcall = squeeze(Pcall(:,1,:));
Pcall = reshape(Pcall,[size(Pcall,1) 4 4]);

for i = 1:4,
    hax(1) = subplot(2,2,i);
    mplot(Lall(:,:,i), Pcall(:,:,i), 'b-|g--|r-|c--','LineWidth',2);
end

    
fexp = catuneven(2,dutydata.fexp);
fexp = reshape(fexp,[size(fexp,1) length(omegarvals) length(zetavals) length(dutyall)]);

figureseries('Floquet exponents vs duty cycle');
lab = {'0.2, 3', '0.2, 0.5', '4, 3', '4, 0.5'};
clf;
for i = 1:2
    hax(i) = subplot(2,1,i);
    plot(dutyall, log(0.5) ./ real(flatten(fexp(i,:,:,:),1:3)), 'o-');
    xlabel('Duty cycle');
    ylabel('t_{1/2} (sec)');
    title(sprintf('Mode %d time constants',i));
    legend(lab, 'Location','NE');
end

dutycyclevals = [0.2 0.36 0.5 0.7];
omegarvals = 2*pi* ([0.5 1 2]);
if (~getvar('nonlindata') || ~inputyn('Use existing nonlinear data?', 'default',true))
    nonlindata = struct([]);
    
    X0 = [L0 - ls0   0   0   0    L0 - ls0  0   0   0   L0   0];
    zeta = zetaold;
    lambda2old = lambda2;
    xmold = xm;
    xpold = xp;
    
    timedWaitBar(0,'Calculating duty cycles');
    a = 1;
    n = length(dutycyclevals) * length(omegarvals) * 3;
    for k = 1:length(dutycyclevals)
        duty = dutycyclevals(k);
        for j = 1:length(omegarvals)
            omegar = omegarvals(j);
            for i = 1:4
                switch i
                    case 1
                        lambda2 = 0;
                        xm = 0;
                        xp = 0;
                    case 2
                        lambda2 = lambda2old;
                        xm = 0;
                        xp = 0;
                    case 3
                        lambda2 = 0;
                        xm = xmold;
                        xp = xpold;
                    case 4
                        lambda2 = lambda2old;
                        xm = xmold;
                        xp = xpold;
                end
                
                fprintf('Duty = %g, OmegaR = %g, Nonlin = %d\n', duty, omegar, i);

                [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, X0, ...
                    'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
                data1.Pc(:,1) = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
                data1.Pc(:,2) = Pc(data1.x(:,5), data1.x(:,6), data1.x(:,8));

                data1 = get_floquet(data1,@(t,x) jfcn(t,x), 150);
                
                data1.dutycycle = duty;
                data1.zeta = zeta;
                data1.omegar = omegar;
                data1.lambda2 = lambda2;
                data1.xm = xm;
                data1.xp = xp;
                
                nonlindata = makestructarray(nonlindata,data1);
                timedWaitBar(a/n);
                a = a+1;
                putvar nonlindata;
            end
        end
    end
end

nonlindata = reshape(nonlindata, [4 length(omegarvals) length(dutycyclevals)]);

Pcall = cat(3, nonlindata.Pc);
Pcall = reshape(Pcall,[size(Pcall,1) 2 4 3 4]);

xx = cat(3, nonlindata.x);
xx = reshape(xx, [size(xx,1) size(xx,2) size(nonlindata)]);
tt = nonlindata(1).t;

fexp = catuneven(2,nonlindata.fexp);
fexp = reshape(fexp,[size(fexp,1) size(nonlindata)]);

figureseries('Floquet exponents vs nonlinearity');
lab = {'flat','FL','FV','FLV'};
clf;
for i = 1:2
    hax(i) = subplot(2,1,i);
    plot(dutycyclevals, log(0.5) ./ real(flatten(fexp(i,:,2,:),1:3)), 'o-');
    xlabel('Duty cycle');
    ylabel('t_{1/2} (sec)');
    title(sprintf('Mode %d time constants',i));
    legend(lab, 'Location','NE');
end

figureseries('Limit cycle vs nonlinearity');
clf;
subplot(2,1,1);
plot(tt,squeeze(Pcall(:,1,:,2,2)));
addplot(tt,squeeze(Pcall(:,2,:,2,2)),'--');

subplot(2,1,2);
plot(tt,squeeze(xx(:,9,:,2,2)));
    
figureseries('Limit cycle vs duty cycle');
clf;
subplot(2,1,1);
plot(tt,squeeze(Pcall(:,1,4,2,:)));

subplot(2,1,2);
plot(tt,squeeze(xx(:,9,4,2,:)));

figureseries('Floquet exponents vs duty cycle2');
clf;
for i = 1:2
    hax(i) = subplot(2,1,i);
    plot(dutycyclevals, log(0.5) ./ real(flatten(fexp(i,4,:,:),1:3)), 'o-');
    xlabel('Duty cycle');
    ylabel('t_{1/2} (sec)');
    title(sprintf('Mode %d time constants',i));
end

    function actval = act(t)
        
        t1 = mod(t,1);
        actval(1,:) = double(t1 < duty);
        actval(2,:) = double(((t1 >= 0.5) & (t1 < 0.5+duty)) | ((t1 >= 0) & (t1 < duty-0.5)));

    end

    function [hx,dhx] = h(x)
        
        exs = exp(x/s);
        hx = s * log(1 + exs);
        if (nargout == 2)
            dhx = exs ./ (1 + exs);
        end
        
    end

    function [l0,dl] = lambdafcn(lc)
        
        l0 = 1 + lambda2 * (lc - lc0).^2;
        if (nargout == 2)
            dl = 2.*lambda2.*(lc - lc0);
            dl(l0 < 0) = 0;
        end
        l0(l0 < 0) = 0;
        
    end

    function [x,dx] = xifcn(vc)
        
        if (nargout == 1)
            x = 1 + xp * h(vc) - xm * h(-vc);
            x(x > ximax) = ximax;
        else
            [hvcp,dhvcp] = h(vc);
            [hvcm,dhvcm] = h(-vc);
            x = 1 + xp * hvcp - xm * hvcm;
            
            dx = xm .* dhvcm + xp .* dhvcp;
            dx(x > ximax) = 0;
            x(x > ximax) = ximax;
        end
        
    end

    function p = Pc(lc, vc, Caf)
        
        p = P0 .* lambdafcn(lc) .* xifcn(vc) .* Caf;
        
    end

    function dx = odefcn(t,x)
        
        lc = x([1 5],:);
        vc = x([2 6],:);
        Ca = x([3 7],:);
        Caf = x([4 8],:);
        L = x(9,:);
        V = x(10,:);
        
        t1 = mod(t,1);
        actval(1,:) = double(t1 < duty);
        actval(2,:) = double(((t1 >= 0.5) & (t1 < 0.5+duty)) | ((t1 >= 0) & (t1 < duty-0.5)));
        
        gact = 1./(1+exp(2-(actval-0.5)/s));
        
        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);
        dlc = vc;

        xi = 1 + xp * h(vc) - xm * h(-vc);
        lambda = 1 + lambda2 * (lc - lc0).^2;
        Pcval = P0 .* lambda .* xi .* Caf;

        dvc(1,:) = 1/mm * (-Pcval(1,:) + mu * (L - lc(1,:) - ls0) - B * vc(1,:));
        dvc(2,:) = 1/mm * (-Pcval(2,:) + mu * (2*L0-L - lc(2,:) - ls0) - B * vc(2,:));
        
        dL = V;
        dV = 1/M * (-mu * (L - lc(1,:) - ls0) + B * vc(1,:) + ...
                     mu * (2*L0-L - lc(2,:) - ls0) - B * vc(2,:) ...
                     - 2 * zeta * omegar * V - omegar^2 * (L - L0));
        
        dx = [dlc(1); dvc(1); dCa(1); dCaf(1); ...
            dlc(2); dvc(2); dCa(2); dCaf(2); ...
            dL; dV];
        
    end

    function J = jfcn(t,x)
         
        lc = reshape(x([1 5]), [1 1 2]);
        vc = reshape(x([2 6]), [1 1 2]);
        Ca = reshape(x([3 7]), [1 1 2]);
        Caf = reshape(x([4 8]), [1 1 2]);
        L = x(9);
        V = x(10);
        
        [lambda,dlambda] = lambdafcn(lc);
        [xi,dxi] = xifcn(vc);
        actval = act(t);
        gact = 1./(1+exp(2-(actval-0.5)/s));
        gact = reshape(gact, [1 1 2]);
        
        J = zeros(10,10);
        
        %each muscle block has an independent 4x4 block in the Jacobian
        %that describes its independent evolution
        zero = zeros(1,1,2);
        one = ones(1,1,2);
        Jmusc = [ ...
            zero, one, zero, zero; ...
            ...
            1/mm .* (-mu - Caf .* P0 .* xi .* dlambda), ...
            1/mm .* (-B - Caf .* P0 .* lambda .* dxi), ...
            zero, ...
            -1/mm .* P0 .* lambda .* xi; ...
            ...
            zero, ...
            zero, ...
            -(1 - Caf) .* k3 - Ca .* k2 .* (1 - gact) + ...
               k2 .*(-Ca - Caf + C - S) .* (1 - gact) - k1 .* gact, ...
            Ca .* k3 + (1 - Caf).*k4 - Caf.*k4 - Ca.*k2 .* (1 - gact) - ...
               k1 .* gact; ...
            ...
            zero, ...
            zero, ...
            (1 - Caf) .* k3, ...
            -Ca .* k3 - (1 - Caf) .* k4 + Caf .* k4
            ];
        
        J(1:4,1:4) = Jmusc(:,:,1);
        J(5:8,5:8) = Jmusc(:,:,2);
        
        
        %spring mass block
        J(9:10,9:10) = [0 1; ...
            1/M*(-2*mu - omegar^2) -2*zeta*omegar/M];
        
        %coupling between muscles and spring mass
        J(10,1:8) = [mu/M  B/M  0  0  -mu/M   -B/M   0   0];
        
        J(1:8,9) = [0; mu/mm; 0; 0; 0; -mu/mm; 0; 0];    
    end

    function stat = odeoutput(t,y, flag)
        if (isempty(flag) && ~isempty(t) && (mod(t(1),1) < 0.0001))
            fprintf('t = %g\n', t(1));
        end
        stat = 0;
    end

end

  

