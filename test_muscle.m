function test_muscle

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
M = 0.05;
B = 10;
xm = 0.4;
xp = 1.33;
xmax = 1.8;
s = 0.1;
L0 = 2.7;
A = 0.125;
phi = 0.1;
T = 1;

act = @(t) mod(t,1) < 0.36;
L = @(t) L0 + A * cos(2*pi/T * (t - phi));

dt = 0.005;
phitest = 0:0.05:0.95;
showphi = [3 7 11 15];

if (~getvar('data') || ~inputyn('Use existing data?', 'default',true))
    data = struct([]);
    for i = 1:length(phitest)
        phi = phitest(i);
        
        L = @(t) L0 + A * cos(2*pi/T * (t - phi));
        X0 = [L(0) - ls0   0   0   0];
        
        [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
            'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
        data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
        
        data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
        data1.L = L;
        data = makestructarray(data,data1);
    end
    putvar data;
end

fxx = cat(4, data.fx);
sgn1 = sign(fxx(1,1,1,:));
fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

fexp = cat(2,data.fexp);

figure(1);
plot(phitest, log(0.5) ./ real(fexp)');
xlabel('Phase');
ylabel('t_{1/2} (sec)');
title('Mode time constants');

t = data(1).t;

figure(2);
clf;
for i = 1:4,
    subplot(2,2,i);
    j = showphi(i);
    plot(data(j).t, fxx(:,:,1,j));
    xlabel('Time (s)');
    title(sprintf('\\phi_{act} = %g',phitest(j)));
end
legend('lc','vc','Ca','Caf','Location','best');

figure(3);
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

figure(4);
clf;
Pcdevall = zeros(length(t),length(phitest),length(phitest));
W0 = zeros(length(phitest),1);
Wdev = zeros(length(phitest),length(phitest));
for i = 1:length(data)
    xbase = data(i).x;
    Pcbase = data(i).Pc;
    
    W0(i) = trapz(-data(i).L(t), Pcbase);
    % plot(-data(i).L(t), Pcbase, 'k-');
    % fprintf('Total work at phase %g = %g\n', phitest(i), W0(i));
    
    dec = exp(data(i).fexp(1)*data(i).t);
    for j = 1:length(phitest),
        a = find(t >= phitest(j),1);
        k = [a:length(t) 1:a-1]';
        dec1 = zeros(size(dec));
        dec1(k) = dec;
        
        dev1 = fxx(:,:,1);
        %dev1 = dev1 / sqrt(sum(dev1(1,:).^2));
        dev1 = bsxfun(@times, dev1, dec1);
        
        xfx1 = xbase + 0.2*dev1;
        Pcdevall(:,i,j) = Pc(xfx1(:,1),xfx1(:,2),xfx1(:,4));
        
        Wdev(i,j) = trapz(-data(i).L(t), Pcdevall(:,i,j));
        %addplot(-data(i).L(t), Pcdevall(:,i,j), 'r-');
        %drawnow;
    end
    %pause;
end

figure(4);
clf;
showphi = [0.1 0.5 0.7];
ishowphi = zeros(size(showphi));
for i = 1:length(showphi)
    j  = find(phitest >= showphi(i), 1);
    ishowphi(i) = j;
end
plot(phitest, W0, phitest, Wdev(:,ishowphi));
xlabel('Activation phase');
ylabel('Work');
legend('steady','perturbed','Location','best');

figure(5);
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

            data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
            data1.L = L;
            L0data = makestructarray(L0data,data1);
        end
    end
    L0data = reshape(L0data,[length(L0test) length(phitest2)]);
    putvar L0data;
end

lcall = cat(3,L0data.x);
lcall = squeeze(lcall(:,1,:));
lcall = reshape(lcall, [size(lcall,1) size(L0data)]);

figure(6);
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

                data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
                data1.L = L;

                NLdata = makestructarray(NLdata,data1);
            end
        end
    end
    NLdata = reshape(NLdata,[length(phitest2) 2 2]);
    putvar NLdata;
end

Pcnl = cat(2,NLdata.Pc);
Pcnl = reshape(Pcnl, [size(Pcnl,1) size(NLdata)]);

figure(7);
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
        
        lc = x(1,:);
        vc = x(2,:);
        Ca = x(3,:);
        Caf = x(4,:);
        
        actval = act(t);
        Lval = L(t);
        
        gact = 1./(1+exp(2-(actval-0.5)/s));
        
        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);
        dlc = vc;
        
        dvc = 1/M * (-Pc(lc,vc,Caf) + mu * (Lval - lc - ls0) - B * vc);
        
        dx = [dlc; dvc; dCa; dCaf];
        
    end

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
        
    end

end

  
