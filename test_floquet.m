function test_floquet(data, ph, dist, varargin)

opt.fixedperiod = false;
opt = parsevarargin(opt, varargin, 4, 'allowunknown');

P = size(data.x,2);

%options for ode
odefields = odeset;
odefields = fieldnames(odefields);
[odeopt,~] = getfieldsonly(opt, odefields);
odeopt = odeset(odeopt);

for j = 1:length(ph)
    t1 = data.t + ph(j)*data.per;
    x0 = deval(data.sol, mod(t1, data.per))';
    for i = 1:length(dist)
        for k = 1:P,
            dec = exp(data.fexp(k)*data.t);
            dxf = fourier_sin_cos(t1, data.fmode(:,:,k), data.per);
            ic = x0(1,:) + dist(i)*dxf(1,:);
            
            if (~opt.fixedperiod)
                odeopt = odeset(odeopt,'Events',@(t,x) istransverse(t,x, fp,data.fcn));
            end
            sol = ode45(data.fcn, [t1(1) t1(end)], ic, odeopt);
            
            xp = deval(sol, t1)';
            dxp = xp - x0;
            
            plot(t1,dxp);
            addplot(t1,dist(i) * bsxfun(@times, dxf, dec), '--');
            pause;
        end
    end
end

function [value,isterminal,direction] = istransverse(t,x, v0,fcn)

ftrn = fcn(t,x);

v = x(1:length(v0));
value = ftrn' * (v-v0);
direction = 1;
isterminal = 1;


