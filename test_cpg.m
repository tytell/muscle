function test_cpg

Gt = 10;
kc = 20;
ta = 0.2;
tx = 0.5;
s = 0.1;

X0 = [1 0 -1 0];
check_jacobian(0,X0', 0.02*ones(10,1), @odefcn, @jfcn);

tinit = [0 15];
odeopt = odeset('RelTol',1e-6); %, 'OutputFcn', @odeoutput);
[t,x] = ode45(@odefcn, tinit, X0, odeopt);

figureseries('Time series');
clf;
hax(1) = subplot(2,1,1);
plot(t,x);
axis tight;

hax(2) = subplot(2,1,2);
plot(t,x(:,[1 3])-x(:,[2 4]));
axis tight;

labellines(hax(1), {'v1','a1','v2','a2'});
labellines(hax(2), {'y1','y2'});

[~,~,data] = get_limit_cycle(@odefcn, 0.005, 1.5, X0, ...
        'Display','iter-detailed', 'fixedperiod',false, 'initialcycles',5, 'TolX',1e-8, 'RelTol',1e-6);
data = get_floquet(data,@(t,x) jfcn(t,x), 100);

figureseries('Limit cycle');
clf;
plot(data.t, data.x);
xlabel('Time');
labellines({'v1','a1','v2','a2'});
axis tight;

figureseries('Floquet modes');
clf;
for i = 1:4,
    subplot(2,2,i);
    plot(data.t, data.fx(:,:,i));
    
    title(sprintf('t_{1/2} = %g', log(0.5) ./ real(data.fexp(i))));
    xlabel('Time');
    labellines({'v1','a1','v2','a2'});
    axis tight;
end

    function dx = odefcn(t,x)
        
        v = x([1 3],:);
        a = x([2 4],:);
        
        hv = s * log(1+exp(v/s));
        
        y = v - a;
        hy = s * log(1+exp(y/s));
        
        dv = 1/tx * (-v + Gt * (1 - v) + ...
            kc * hy([2 1],:) .* (-1 - v));
        
        da = 1/ta * (hv - a);
        
        dx = [dv(1); da(1); dv(2); da(2)];
        
    end

    function J = jfcn(t,x)
        
        v = reshape(x([1 3]), [1 1 2]);
        evs = exp(v/s);
        dhv = evs ./ (1 + evs);
        
        a = reshape(x([2 4]), [1 1 2]);
        y = v - a;

        eys = exp(y/s);
        hy = s * log(1 + eys);
        dhy = eys ./ (1 + eys);
             
        %block diagonal 2x2
        zero = zeros(1,1,2);
        one = ones(1,1,2);
        Jd = [ ...
            1/tx * (-1 - Gt - kc * hy(1,1,[2 1])),   zero; ...
            dhv / ta,                                -1/ta * one];
        %off diagonal 2x2
        Jo = [ ...
            1/tx * kc * (-1 - v) .* dhy(1,1,[2 1]),   -1/tx * kc * (-1 - v) .* dhy(1,1,[2 1]); ...
            zero,                                     zero];
        
        J = [Jd(:,:,1) Jo(:,:,1); ...
            Jo(:,:,2) Jd(:,:,2)];
            
    end
            
end