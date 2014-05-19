function test_muscle_mass

%optimized parameters:
%2: m = 0.0542; b = 0.2802; lc0 = 0.9678; k1 = 6.7281; k2 = 23.2794; k30 = 51.3537; k40 = 19.3801; km1 = 17.5804; km2 = 6.0156    ->> sum(dx^2) = 6.056118

filename = 'test_muscle_mass.mat';
quiet = true;
doplot = false;

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
    

  
