function test_two_muscles3

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
xmax = 1.8;
s = 0.1;
dphi = 0.5;
T = 1;
duty = 0.36;

L0 = 2.7;
M = 0.5;
zeta = 0.5;
omegar = 2;

dt = 0.005;

omegarvals = 2*pi* ([0.3:0.05:1 1.2 1.5 2]);
showfreq = [1 7 15 18];

X0 = [L0 - ls0   0   0   0    L0 - ls0  0   0   0   L0   0];
check_jacobian(0,X0', 0.02*ones(10,1), @odefcn, @jfcn);

tinit = [0 15*T];
odeopt = odeset('RelTol',1e-6); %, 'OutputFcn', @odeoutput);
[t,x] = ode45(@odefcn, tinit, X0, odeopt);

figureseries('Time series');
clf;
hax(1) = subplot(3,1,1);
plot(t,x(:,1:4));
ylabel('Left muscle');
axis tight;

hax(2) = subplot(3,1,2);
plot(t,x(:,5:8),'--');
ylabel('Right muscle');
axis tight;

hax(3) = subplot(3,1,3);
addmplot(t,x(:,9:10),'k-|k:', 'LineWidth',2);
axis tight;

linkaxes(hax, 'x');
set(hax, 'XLim',[10 15]);

labellines(hax(1), {'lc','vc','Ca','Caf'}, 'location',[13.25 13.25 14.4 14.5]);
labellines(hax(3), {'L','V'});

if (~getvar('freqdata') || ~inputyn('Use existing data?', 'default',true))
    freqdata = struct([]);
    
    for i = 1:length(omegarvals)
        omegar = omegarvals(i);
        
        [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',15, 'TolX',1e-8, 'RelTol',1e-6);
        data1.Pc(:,1) = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
        data1.Pc(:,2) = Pc(data1.x(:,5), data1.x(:,6), data1.x(:,8));
        data1.omegar = omegar;
        data1.zeta = zeta;
        data1 = get_floquet(data1,@(t,x) jfcn(t,x), 150);
        freqdata = makestructarray(freqdata,data1);
    end
    
    putvar freqdata;
end

figureseries('Effect of resonant frequency');
xx = cat(3, freqdata.x);
Pcall = cat(3, freqdata.Pc);

subplot(1,2,1);
imagesc(freqdata(1).t, 1-omegarvals/(2*pi), squeeze(xx(:,9,:))');
ylabel('Frequency difference f_{act} - f_{res} (Hz)');
xlabel('Time (sec)');
hcol = colorbar;
ylabel(hcol, 'L');

subplot(1,2,2);
imagesc(freqdata(1).t, 1-omegarvals/(2*pi), squeeze(Pcall(:,1,:))' / P0);
ylabel('Frequency difference f_{act} - f_{res} (Hz)');
xlabel('Time (sec)');
hcol = colorbar;
ylabel(hcol, 'Force ratio P_{c,1} / P_0');

fxx = cat(4, freqdata.fx);
sgn1 = sign(fxx(1,1,1,:));
fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

fexp = cat(2,freqdata.fexp);

figureseries('Floquet exponents vs resonant frequency');
clf;
plot(1-omegarvals/(2*pi), log(0.5) ./ real(fexp(1,:))');
xlabel('Frequency difference f_{act} - f_{res} (Hz)');
ylabel('t_{1/2} (sec)');
title('Mode one time constants');

figureseries('Floquet modes vs resonant frequency');
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
%legend('lc','vc','Ca','Caf','L','V','Location','best');

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

figureseries('Damping and resonance');
clf;
for i = 1:2
    hax(i) = subplot(2,1,i);
    plot(1-omegarvals/(2*pi), log(0.5) ./ real(squeeze(fexp(i,:,:))));
    xlabel('Frequency difference f_{act} - f_{res} (Hz)');
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
    
    title(sprintf('f_{act} - f_{res} = %g, \\zeta = %g, t_{1/2} = %g', ...
        1-omegarvals(i)/(2*pi), zetavalall(j), log(0.5) ./ real(fexp(1,i,j))));
end
labellines(hln, {'lc','vc','Ca','Caf','L','V'});

zetaold = zetavals;
omegarold = omegarvals;
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
                data1.dutycycle = dutycycle;
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

isomega = ismember(omegarold,omegarvals);
iszeta = ismember(zetaold,zetavals);
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
clf;
for i = 1:2
    hax(i) = subplot(2,1,i);
    plot(dutyall, log(0.5) ./ real(flatten(fexp(i,:,:,:),1:3)));
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
        end
        
    end

    function [x,dx] = xifcn(vc)
        
        if (nargout == 1)
            x = 1 + xp * h(vc) - xm * h(-vc);
        else
            [hvcp,dhvcp] = h(vc);
            [hvcm,dhvcm] = h(-vc);
            x = 1 + xp * hvcp - xm * hvcm;
            
            dx = xm .* dhvcm + xp .* dhvcp;
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

  

