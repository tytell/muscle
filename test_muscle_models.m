function test_muscle_models

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

X0 = [L0 - ls0;   0;   0];
Xp0 = [0;  0;  0];

tinit = [0 2*T];
odeopt = odeset('RelTol',1e-6, 'MaxStep',0.1, 'Jacobian',@fjac);
[t,x] = ode15i(@odefcn, tinit, X0,Xp0, odeopt);

figureseries('Time series');
clf;
plot(t,x);

[~,~,data] = get_limit_cycle(@odefcn, 0.005, T, [X0 Xp0]', 'odetype','implicit', ...
    'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6, ...
    'MaxStep',0.1, 'Jacobian',@fjac);
data.Pc = Pcfcn(data.x(:,1), data.xp(:,1), data.x(:,3));

data = get_floquet(data,@fjac, 100, 'odetype','implicit');

plot(data.t, data.fx(:,:,1));

test_floquet(data, 0:0.25:0.75, [0.2 0.1 0.05], 'fixedperiod',true, ...
    'Jacobian',@fjac, 'MaxStep',0.1, 'odetype','implicit');

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
                
        t1 = mod(t,1);
        actval = double((t1 >= 0.1) & (t1 < 0.1+duty)); %(sin(pi * (t/T-dphi))).^2;        
        gact = actval; %1./(1+exp(-actval/s));
        
        L = L0 + A*sin(2*pi * (t/T - dphi));

        dCaf = (k3 * Ca - k4 * Caf) .* (1 - Caf);
        dCa = (k4 * Caf - k3 * Ca) .* (1 - Caf) + ...
            gact .* k1 .* (C - Ca - Caf) + ...
            (1 - gact) .* k2 .* Ca .* (C - S - Ca - Caf);

        xi = 1 + xip * h(vc) - xim * h(-vc);
        lambda = 1 + lambda2 * (lc - lc0).^2;
        Pcval = P0 .* lambda .* xi .* Caf;

        lceqn = -B*vc + mu*(L - lc - ls0) - Pcval;
        
        eqn = [lceqn; Cap-dCa; Cafp-dCaf];
        
    end

    function [dfdx,dfdxp] = fjac(t,x,xp)
        
        lc = x(1);
        Ca = x(2);
        Caf = x(3);
        
        vc = xp(1);

        [lambda,dlambda] = lambdafcn(lc);
        [xi,dxi] = xifcn(vc);

        t1 = mod(t,1);
        actval = double((t1 >= 0.1) & (t1 < 0.1+duty)); %(sin(pi * (t/T-dphi))).^2;        
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
