function test_muscle_mass

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
A = 0.125;
phi = 0.1;
T = 1;

L0 = 2.7;
M = 0.5;
zeta = 0.5;
omegar = 2;

act = @(t) mod(t,1) < 0.36;

dt = 0.005;

omegarvals = 2*pi* ([0.3:0.05:1 1.2 1.5 2]);
showfreq = [1 7 15 18];

if (~getvar('freqdata') || ~inputyn('Use existing data?', 'default',true))
    freqdata = struct([]);
    
    for i = 1:length(omegarvals)
        omegar = omegarvals(i);
        
        X0 = [L0 - ls0   0   0   0    L0   0];

        [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',15, 'TolX',1e-8, 'RelTol',1e-6);
        data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
        data1.omegar = omegar;
        data1.zeta = zeta;
        data1 = get_floquet(data1,@(t,x) jfcn(t,x), 150, 'neigenvalues',8);
        freqdata = makestructarray(freqdata,data1);
    end
    
    putvar freqdata;
end

figure(1);
xx = cat(3, freqdata.x);
Pcall = cat(2, freqdata.Pc);

subplot(1,2,1);
imagesc(freqdata(1).t, 1-omegarvals/(2*pi), squeeze(xx(:,5,:))');
ylabel('Frequency difference f_{act} - f_{res} (Hz)');
xlabel('Time (sec)');
hcol = colorbar;
ylabel(hcol, 'L');

subplot(1,2,2);
imagesc(freqdata(1).t, 1-omegarvals/(2*pi), Pcall' / P0);
ylabel('Frequency difference f_{act} - f_{res} (Hz)');
xlabel('Time (sec)');
hcol = colorbar;
ylabel(hcol, 'Force ratio P_c / P_0');

fxx = cat(4, freqdata.fx);
sgn1 = sign(fxx(1,1,1,:));
fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

fexp = cat(2,freqdata.fexp);

figure(2);
clf;
plot(1-omegarvals/(2*pi), log(0.5) ./ real(fexp(1,:))');
xlabel('Frequency difference f_{act} - f_{res} (Hz)');
ylabel('t_{1/2} (sec)');
title('Mode one time constants');

figure(3);
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

zetaold = zeta;
omegarold = omegarvals;

zetavals = [0.2 1 2 4];
omegarvals = 2*pi* ([0.3 0.5 0.8 1 1.2 1.5 2]);

if (~getvar('dampdata') || ~inputyn('Use existing data?', 'default',true))
    dampdata = struct([]);
    
    X0 = [L0 - ls0   0   0   0    L0   0];
    
    timedWaitBar(0,'Calculating damping and resonance');
    a = 1;
    n = length(omegarvals) * length(zetavals);
    for j = 1:length(zetavals)
        zeta = zetavals(j);
        for i = 1:length(omegarvals)
            omegar = omegarvals(i);
            
            [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, X0, ...
                'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
            data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
            
            data1 = get_floquet(data1,@(t,x) jfcn(t,x), 150, 'neigenvalues',8);
            data1.zeta = zeta;
            data1.omegar = omegar;
            dampdata = makestructarray(dampdata,data1);
            timedWaitBar(a/n,'Calculating damping and resonance');
            a = a+1;
        end
    end
    putvar dampdata;
end

isomega = ismember(omegarold,omegarvals);
islowzeta = zetavals < zetaold;

dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)]);
dampdata = makestructarray(dampdata(:,islowzeta), freqdata(isomega), dampdata(:,~islowzeta));
dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)+1]);
fexp = cat(2,dampdata.fexp);
fexp = reshape(fexp,[size(fexp,1) length(omegarvals) length(zetavals)+1]);

figure(4);
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

% 
% figure(3);
% clf;
% for i = 1:4,
%     subplot(2,2,i);
%     j = showphi(i);
%     h1 = plot(data(j).t, data(j).x, 'k--');
%     h2 = addplot(data(j).t, data(j).x + 0.1*fxx(:,:,1,j));
%     xlabel('Time (s)');
%     title(sprintf('\\phi_{act} = %g',phitest(j)));
% end
% legend([h1(1); h2], 'steady','lc','vc','Ca','Caf','Location','best');
% 
% figure(4);
% clf;
% Pcdevall = zeros(length(t),length(phitest),length(phitest));
% W0 = zeros(length(phitest),1);
% Wdev = zeros(length(phitest),length(phitest));
% for i = 1:length(data)
%     xbase = data(i).x;
%     Pcbase = data(i).Pc;
%     
%     %W0(i) = trapz(-data(i).L(t), Pcbase);
%     % plot(-data(i).L(t), Pcbase, 'k-');
%     % fprintf('Total work at phase %g = %g\n', phitest(i), W0(i));
%     
%     dec = exp(data(i).fexp(1)*data(i).t);
%     for j = 1:length(phitest),
%         a = find(t >= phitest(j),1);
%         k = [a:length(t) 1:a-1]';
%         dec1 = zeros(size(dec));
%         dec1(k) = dec;
%         
%         dev1 = fxx(:,:,1);
%         %dev1 = dev1 / sqrt(sum(dev1(1,:).^2));
%         dev1 = bsxfun(@times, dev1, dec1);
%         
%         xfx1 = xbase + 0.2*dev1;
%         Pcdevall(:,i,j) = Pc(xfx1(:,1),xfx1(:,2),xfx1(:,4));
%         
%         %Wdev(i,j) = trapz(-data(i).L(t), Pcdevall(:,i,j));
%         %addplot(-data(i).L(t), Pcdevall(:,i,j), 'r-');
%         %drawnow;
%     end
%     %pause;
% end
% 
% figure(4);
% clf;
% showphi = [0.1 0.5 0.7];
% ishowphi = zeros(size(showphi));
% for i = 1:length(showphi)
%     j  = find(phitest >= showphi(i), 1);
%     ishowphi(i) = j;
% end
% plot(phitest, W0, phitest, Wdev(:,ishowphi));
% xlabel('Activation phase');
% ylabel('Work');
% legend('steady','perturbed','Location','best');
% 
% figure(5);
% clf;
% contourf(phitest, phitest, bsxfun(@minus, Wdev, W0)' / max(abs(W0)));
% hcol = colorbar;
% 
% xlabel('Phase of activation');
% ylabel('Phase of perturbation');
% 
% ylabel(hcol, 'Fractional change');
% 
% 
% L0test = lc0 + ls0 + [-2*A 0 2*A];
% phitest2 = 0:0.2:0.8;
% if (~getvar('L0data') || ~inputyn('Use existing data?', 'default',true))
%     L0data = struct([]);
%     for i = 1:length(phitest2)
%         phi = phitest2(i);
%         for j = 1:length(L0test)
%             L0 = L0test(j);
% 
%             % L = @(t) L0 + A * cos(2*pi/T * (t - phi));
%             X0 = [L0 - ls0   0   0   0   L0  0];
% 
%             [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
%                 'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
%             data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
% 
%             data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
%             L0data = makestructarray(L0data,data1);
%         end
%     end
%     L0data = reshape(L0data,[length(L0test) length(phitest2)]);
%     putvar L0data;
%     
%     %reset L0
%     L0 = 2.7;
% end
% 
% lcall = cat(3,L0data.x);
% lcall = squeeze(lcall(:,1,:));
% lcall = reshape(lcall, [size(lcall,1) size(L0data)]);
% 
% figure(6);
% clf;
% subplot(2,2,1);
% lc1 = 0.7*lc0:0.01:1.3*lc0;
% plot(lc1, lambdafcn(lc1), 'k--');
% addplot(squeeze(lcall(:,1,:)), squeeze(lambdafcn(lcall(:,1,:))),'b-', ...
%     squeeze(lcall(:,2,:)), squeeze(lambdafcn(lcall(:,2,:))),'g-', ...
%     squeeze(lcall(:,3,:)), squeeze(lambdafcn(lcall(:,3,:))),'r-', ...
%     'LineWidth',2);
% axis tight;
% xlabel('lc');
% ylabel('\lambda');
% 
% subplot(2,2,2);
% plot(t, squeeze(lcall(:,1,:)),'b-', ...
%     t, squeeze(lcall(:,2,:)),'g-', ...
%     t, squeeze(lcall(:,3,:)),'r-');
% xlabel('Time (s)');
% ylabel('lc');
% 
% subplot(2,1,2);
% fx = cat(3,L0data.fexp);
% fx = reshape(fx,[4 size(L0data)]);
% plot(phitest2,squeeze(fx(1,:,:)),'o-');
% xlabel('Activation phase');
% ylabel('Mode 1 exponent');
% 
% if (~getvar('NLdata') || ~inputyn('Use existing data?', 'default',true))
%     NLdata = struct([]);
%     for i = 1:2
%         switch i
%             case 1
%                 lambda2 = -2.23;
%             case 2
%                 lambda2 = 0;
%         end
%         for j = 1:2
%             switch j
%                 case 1
%                     xm = 0.4;
%                     xp = 1.33;
%                 case 2
%                     xm = 0;
%                     xp = 0;
%             end
% 
%             for k = 1:length(phitest2)
%                 phi = phitest2(k);
%                 
%                 %L = @(t) L0 + A * cos(2*pi/T * (t - phi));
%                 X0 = [L0 - ls0   0   0   0   L0  0];
% 
%                 [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
%                     'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
%                 data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
% 
%                 data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
% 
%                 NLdata = makestructarray(NLdata,data1);
%             end
%         end
%     end
%     NLdata = reshape(NLdata,[length(phitest2) 2 2]);
%     putvar NLdata;
% 
% 
%     lambda2 = -2.23;
%     xm = 0.4;
%     xp = 1.33;
% end
% 
% Pcnl = cat(2,NLdata.Pc);
% Pcnl = reshape(Pcnl, [size(Pcnl,1) size(NLdata)]);
% 
% figure(7);
% clf;
% for i = 1:5,
%     subplot(2,3,i);
%     plot(t, flatten(Pcnl(:,i,:,:),2:4));
%     xlabel('Time (s)');
%     ylabel('Force (mN)');
%     
%     if (i == 3)
%         legend('normal','fv=0','fl=0','both', 'location','best');
%     end
% end
% 
% 
% subplot(2,3,6);
% fx = cat(3,NLdata.fexp);
% fx = reshape(fx,[4 size(NLdata)]);
% plot(phitest2,squeeze(fx(1,:,:)),'o-');
% xlabel('Activation phase');
% ylabel('Mode 1 exponent');
% 
% Btest = [5 10 15 20];
% if (~getvar('Bdata') || ~inputyn('Use existing B data?', 'default',true))
%     Bdata = struct([]);
%     for i = 1:length(phitest2)
%         phi = phitest2(i);
%         for j = 1:length(Btest)
%             B = Btest(j);
%             
%             %L = @(t) L0 + A * cos(2*pi/T * (t - phi));
%             X0 = [L0 - ls0   0   0   0   L0   0];
% 
%             [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
%                 'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
%             data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
% 
%             data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
%             data1.B = B;
%             Bdata = makestructarray(Bdata,data1);
%         end
%     end
%     Bdata = reshape(Bdata,[length(Btest) length(phitest2)]);
%     putvar Bdata;
%     
%     B = 10;
% end
% 
% figure(8);
% fx = cat(3,Bdata.fexp);
% fx = reshape(fx,[4 size(Bdata)]);
% plot(Btest,log(0.5) ./ real(squeeze(fx(1,:,:))),'o-');
% 
% xlabel('Damping coefficient');
% ylabel('t_{1/2} (sec)');
% title('Mode one time constant vs damping');
% 
% k3test = 0.6:0.1:1.4;
% k4test = [0.8 1 1.2];
% if (~getvar('k34data') || ~inputyn('Use existing k3 k4 data?', 'default',true))
%     phi = 0;
%     k30 = k3;
%     k40 = k4;
%     
%     %L = @(t) L0 + A * cos(2*pi/T * (t - phi));
%     X0 = [L0 - ls0   0   0   0   L0   0];
%     
%     k34data = struct([]);
%     for i = 1:length(k3test)
%         k3 = k30*k3test(i);
%         for j = 1:length(k4test)
%             k4 = k40*k4test(j);
%             
%             [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
%                 'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
%             data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
% 
%             data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
%             data1.k3 = k3;
%             data1.k4 = k4;
%             k34data = makestructarray(k34data,data1);
%         end
%     end
%     k34data = reshape(k34data,[length(k4test) length(k3test)]);
%     putvar k34data;
%     
%     k3 = k30;
%     k4 = k40;
% end
% 
% figure(9);
% fx = cat(3,k34data.fexp);
% fx = reshape(fx,[4 size(k34data)]);
% plot(k3test*k3,log(0.5) ./ real(squeeze(fx(1,:,:))),'o-');
% 
% lab = cell(size(k4test));
% for i = 1:length(k4test)
%     lab{i} = sprintf('%.2g',k4test(i)*k4);
%     if (i == 1)
%         lab{i} = ['k4 = ' lab{i}];
%     end
% end
% labellines(lab ,'rotation',0);
% 
% xlabel('k3');
% ylabel('t_{1/2} (sec)');
% title('Mode one time constant vs k3 and k4');
% 

    function [hx,dhx] = h(x)
        
        exs = exp(x/s);
        hx = s * log(1 + exs);
        if (nargout == 2)
            dhx = exs ./ (1 + exs);
        end
        
    end

    function [hl,dl] = lambdafcn(lc)
        
        l0 = 1 + lambda2 * (lc - lc0).^2;
        if (nargout == 1)
            hl = h(l0);
        else
            [hl,dhl] = h(l0);
            dl = 2.*lambda2.*(lc - lc0) .* dhl;
        end
        
    end

    function [x,dx] = xifcn(vc)
        
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
        
        p = P0 .* lambdafcn(lc) .* xifcn(vc) .* Caf;
        
    end

    function dx = odefcn(t,x)
        
        lc = x(1,:);
        vc = x(2,:);
        Ca = x(3,:);
        Caf = x(4,:);
        L = x(5,:);
        V = x(6,:);
        
        actval = act(t);
        %Lval = L(t);
        
        gact = 1./(1+exp(2-(actval-0.5)/s));
        
        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);
        dlc = vc;
        
        dvc = 1/mm * (-Pc(lc,vc,Caf) + mu * (L - lc - ls0) - B * vc);
        
        dL = V;
        dV = 1/M * (-mu * (L - lc - ls0) + B * vc - 2 * zeta * omegar * V - omegar^2 * (L - L0));
        
        dx = [dlc; dvc; dCa; dCaf; dL; dV];
        
    end

    function J = jfcn(t,x)
        
        % From Mathematica:
        % [0,1,0,0,0,0;M.^(-1).*((-1).*\[Mu]+(-1).*Caf.*P0.*xi(vc).*d( ...
        % lambda)(lc)),M.^(-1).*((-1).*B+(-1).*Caf.*P0.*lambda(lc).*d( ...
        % xi)(vc)),0,(-1).*M.^(-1).*P0.*lambda(lc).*xi(vc),M.^(-1).*\[Mu], ...
        % 0;0,0,(-1).*(1+(-1).*Caf).*k3+(-1).*Ca.*k2.*(1+(-1).*g(act( ...
        % t)))+k2.*((-1).*Ca+(-1).*Caf+Cb+(-1).*Sb).*(1+(-1).*g(act(t) ...
        % ))+(-1).*k1.*g(act(t)),Ca.*k3+(1+(-1).*Caf).*k4+(-1).*Caf.* ...
        % k4+(-1).*Ca.*k2.*(1+(-1).*g(act(t)))+(-1).*k1.*g(act(t)),0, ...
        % 0;0,0,(1+(-1).*Caf).*k3,(-1).*Ca.*k3+(-1).*(1+(-1).*Caf).* ...
        % k4+Caf.*k4,0,0;0,0,0,0,0,1;0,0,0,0,(-1).*km,(-1).*bm];
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
        L = x(5,:);
        V = x(6,:);
        
        [l,dl] = lambdafcn(lc);
        [xi,dxi] = xifcn(vc);
        actval = act(t);
        gact = 1./(1+exp(2-(actval-0.5)/s));
        
        J = [0,1,0,0,0,0; ...
            ...
            1/mm .* (-mu - Caf .* P0 .* xi .* dl), ...
            1/mm .* (-B - Caf .* P0 .* l .* dxi), ...
            0, ...
            -1/mm .* P0 .* l .* xi, ...
            mu/mm, 0; ...
            ...
            0, 0, ...
            -(1 - Caf) .* k3 - Ca .* k2 .* (1 - gact) + ...
               k2 .*(-Ca - Caf + C - S) .* (1 - gact) - k1 .* gact, ...
            Ca .* k3 + (1 - Caf).*k4 - Caf.*k4 - Ca.*k2 .* (1 - gact) - ...
               k1 .* gact, ...
            0, 0;...
            ...
            0, 0, (1 - Caf).*k3, -Ca.*k3 - (1 - Caf).*k4 + Caf.*k4, 0, 0;
            0, 0, 0, 0, 0, 1; ...
            mu/M, B/M, 0, 0, 1/M*(-mu-omegar^2), -2*zeta*omegar/M];
        
    end

end

  
