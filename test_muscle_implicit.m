function test_muscle_implicit

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
B = 100;
xim = 0.4;
xip = 1.33;
s = 0.1;
dphi = 0;
duty = 0.36;

L0 = 2.7;
A = 0.05*L0;
T = 1;

L = @(t) L0 + A*sin(2*pi * (t/T - dphi));
act = @(t) double(mod(t,1) < duty); 

X0 = [L0 - ls0   0   0];
Xp0 = get_ic_deriv(@odefcn, @fjac, X0);

tinit = [0 2*T];
odeopt = odeset('RelTol',1e-6, 'MaxStep',0.1, 'Jacobian',@fjac);
[t,x] = ode15i(@odefcn, tinit, X0',Xp0', odeopt);

figureseries('Time series');
clf;
plot(t,x);

phitest = 0:0.05:0.95;
showphi = [1 6 11 17];

if (~getvar('data') || ~inputyn('Use existing data?', 'default',true))
    data = struct([]);
    for i = 1:length(phitest)
        L = @(t) L0 + A*sin(2*pi * (t/T - phitest(i)));
        X0 = [L(0) - ls0   0   0];
        Xp0 = get_ic_deriv(@odefcn, @fjac, X0);
        
        [~,~,data1] = get_limit_cycle(@odefcn, 0.005, T, [X0; Xp0], ...
            'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, ...
            'RelTol',1e-6, 'MaxStep',0.1, 'Jacobian',@fjac, 'odetype','implicit');
        data1.Pc = Pcfcn(data1.x(:,1), data1.xp(:,1), data1.x(:,3));
        
        data1 = get_floquet(data1,@fjac, 100, 'odetype','implicit');
        data1.L = L;
        data1.dphi = phitest(i);
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
dt = data(1).t(2) - data(1).t(1);

figureseries('Floquet modes vs. phi');
clf;
for i = 1:4,
    subplot(2,2,i);
    j = showphi(i);
    plot(data(j).t, fxx(:,:,1,j));
    xlabel('Time (s)');
    title(sprintf('\\phi_{act} = %g',phitest(j)));
end
legend('lc','Ca','Caf','Location','best');

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
legend([h1(1); h2], 'steady','lc','Ca','Caf','Location','best');

Pcdevall = zeros(length(t),length(phitest),length(phitest),2);
W0 = zeros(length(phitest),1);
Wdev = zeros(length(phitest),length(phitest),2);
for f = 1:2,
    for i = 1:length(data)
        xbase = data(i).x;
        xpbase = data(i).xp;
        Pcbase = data(i).Pc;

        W0(i) = trapz(-data(i).L(t), Pcbase);
        % plot(-data(i).L(t), Pcbase, 'k-');
        % fprintf('Total work at phase %g = %g\n', phitest(i), W0(i));

        dec = exp(data(i).fexp(f)*data(i).t);
        ddec = data(i).fexp(f) * exp(data(i).fexp(f)*data(i).t);

        fx1 = fxx(:,:,f);
        dfx1 = deriv(dt, fx1([end-1 1:end 2],:));
        dfx1 = dfx1(2:end-1,:);
        for j = 1:length(phitest),
            a = find(t >= phitest(j),1);
            k = [a:length(t) 1:a-1]';
            dec1 = zeros(size(dec));
            dec1(k) = dec;
            ddec1 = zeros(size(dec));
            ddec1(k) = ddec;

            %dev1 = dev1 / sqrt(sum(dev1(1,:).^2));
            dev1 = bsxfun(@times, fx1, dec1);
            ddev1 = bsxfun(@times, fx1, ddec1) + bsxfun(@times, dfx1, dec1);
            
            xfx1 = xbase + 0.2*dev1;
            xfx1p = xpbase + 0.2*ddev1;
            Pcdevall(:,i,j,f) = Pcfcn(xfx1(:,1),xfx1p(:,1),xfx1(:,3));

            Wdev(i,j,f) = trapz(-data(i).L(t), Pcdevall(:,i,j,f));
            %addplot(-data(i).L(t), Pcdevall(:,i,j), 'r-');
            %drawnow;
        end
        %pause;
    end
end

figureseries('Effect of perturbations');
clf;
for f = 1:2
    subplot(2,1,f);
    showphi = 1;
    plot(data(showphi).t, data(showphi).Pc, 'k--','LineWidth',2);
    for i = 1:4:length(phitest),
        addplot(data(showphi).t, Pcdevall(:,showphi,i,f), 'r-');
    end
end

figureseries('Change in work');
clf;
showphi = 1:5:length(phitest);
for f = 1:2
    subplot(2,1,f);
    
    plot(phitest, bsxfun(@minus,Wdev(showphi,:,f),W0(showphi)) / max(abs(W0)));
    xlabel('Perturbation phase');
    ylabel('Work');

    lab = cell(size(showphi));
    for i = 1:length(showphi)
        lab{i} = num2str(phitest(showphi(i)));
    end

    labellines(lab, 'location',[0.6 0.6 0.65 0.45],'rotation',0);
end

figureseries('Change in work contour');
clf;
contourf(phitest, phitest, bsxfun(@minus, Wdev(:,:,1), W0)' / max(abs(W0)));
hcol = colorbar;

xlabel('Phase of activation');
ylabel('Phase of perturbation');



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
            x = 1 + xip * h(vc) - xim * h(-vc);
        else
            [hvcp,dhvcp] = h(vc);
            [hvcm,dhvcm] = h(-vc);
            x = 1 + xip * hvcp - xim * hvcm;
            
            dx = xim .* dhvcm + xip .* dhvcp;
        end
        
    end

    function Pc = Pcfcn(lc, vc, Caf)
        
        Pc = P0 .* lambdafcn(lc) .* xifcn(vc) .* Caf;
        
    end

    function eqn = odefcn(t,x,xp)
        
        lc = x(1,:);
        Ca = x(2,:);
        Caf = x(3,:);
        
        vc = xp(1,:);
        Cap = xp(2,:);
        Cafp = xp(3,:);

        Lval = L(t);
        
        actval = act(t);
        gact = actval; %1./(1+exp(-actval/s));
        
        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);

        xi = 1 + xip * h(vc) - xim * h(-vc);
        lambda = 1 + lambda2 * (lc - lc0).^2;
        Pcval = P0 .* lambda .* xi .* Caf;

        lceqn = -B*vc + mu*(Lval - lc - ls0) - Pcval;
        
        eqn = [lceqn; Cap-dCa; Cafp-dCaf];
        
    end

    function [dfdx,dfdxp] = fjac(t,x,xp)
        
        lc = x(1);
        Ca = x(2);
        Caf = x(3);
        
        vc = xp(1);

        [lambda,dlambda] = lambdafcn(lc);
        [xi,dxi] = xifcn(vc);

        actval = act(t);
        gact = actval; %1./(1+exp(-actval/s));

        dfdx = [-mu - Caf*P0*xi*dlambda, 0, -P0*lambda*xi; ...
            0, ...
              (1-Caf)*k3 + Ca*k2*(1 - gact) - k2*(C - S - Ca - Caf)*(1 - gact) + k1*gact, ...
              -Ca*k3 - (1 - Caf)*k4 + Caf*k4 + Ca*k2*(1 - gact) + k1*gact; ...
            0, -(1-Caf)*k3, Ca*k3 + (1-Caf)*k4 - Caf*k4];
        
        dfdxp = [-B - Caf*P0*lambda*dxi, 0, 0; ...
            0, 1, 0; ...
            0, 0, 1];
    end
end
