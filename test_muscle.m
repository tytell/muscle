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
t1 = 0:dt:2*T;

phitest = [0.1 0.3 0.5 0.7];
data = struct([]);
for i = 1:length(phitest)
    phi = phitest(i);

    L = @(t) L0 + A * cos(2*pi/T * (t - phi));
    X0 = [L(0) - ls0   0   0   0];

    [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x), 0.005, T, X0, ...
        'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);
    data1.Pc = Pc(data1.x(:,1), data1.x(:,2), data1.x(:,4));
    
    data1 = get_floquet(data1,@(t,x) jfcn(t,x), 100, 'neigenvalues',4);
    data = makestructarray(data,data1);
end

fexp = cat(2,data.fexp);
plot(phitest, fexp');

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

  
