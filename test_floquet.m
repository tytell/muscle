function test_floquet(data, ph, dist, varargin)

opt.fixedperiod = false;
opt.odetype = 'direct';
opt = parsevarargin(opt, varargin, 4, 'allowunknown');

P = size(data.x,2);

%options for ode
odefields = odeset;
odefields = fieldnames(odefields);
[odeopt,~] = getfieldsonly(opt, odefields);
odeopt = odeset(odeopt);

if strcmp(opt.odetype, 'implicit')
    [~,dfdxp] = opt.Jacobian(0,data.x(1,:)', ones(P,1));
    isimp = any(dfdxp - eye(P) ~= 0);
end

for j = 1:length(ph)
    t1 = data.t + ph(j)*data.per;
    [x0,xp0] = deval(data.sol, mod(t1, data.per));
    x0 = x0';
    xp0 = xp0';
    for i = 1:length(dist)
        for k = 1:P,
            if ~isfinite(data.fexp(k))
                continue
            end
            dec = exp(data.fexp(k)*data.t);
            dxf = fourier_sin_cos(t1, data.fmode(:,:,k), data.per);
            ic = x0(1,:) + dist(i)*dxf(1,:);
            
            if (~opt.fixedperiod)
                odeopt = odeset(odeopt,'Events',@(t,x) istransverse(t,x, fp,xp0(1,:)));
            end
            switch opt.odetype
                case 'direct'
                    sol = ode45(data.fcn, [t1(1) t1(end)], ic, odeopt);
                    
                case 'implicit'
                    icd = zeros(size(ic));
                    eqn = data.fcn(0,ic',icd');
                    icd(~isimp) = -eqn(~isimp);
                    [ic,icd] = decic(data.fcn, 0, ic,zeros(size(ic)), icd,~isimp);        

                    if (~opt.fixedperiod)
                        odeopt = odeset(odeopt, 'Events',@(t,x,xp) istransverse(t,x, fp,xp0(1,:)));
                    end
                    sol = ode15i(data.fcn, t1([1 end]), ic, icd, odeopt);
            end
            
            [xd, xdp] = deval(sol, t1);
            xd = xd';
            xdp = xdp';
            dxd = xd - x0;
            
            plot(t1,dxd);
            addplot(t1,dist(i) * bsxfun(@times, dxf, dec), '--');
            pause;
        end
    end
end

function [value,isterminal,direction] = istransverse(t,x, v0,ftrn)

v = x(1:length(v0));
value = ftrn' * (v-v0);
direction = 1;
isterminal = 1;


