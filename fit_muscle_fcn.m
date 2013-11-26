function [dx,Pcmb,data] = fit_muscle_fcn(fitpar, odefcn, par)

par.mm = fitpar(1);
par.b = fitpar(2);
par.lc0 = fitpar(3);
par.k1 = fitpar(4);  
par.k2 = fitpar(5); 
par.k30 = fitpar(6);
par.k40 = fitpar(7); 
par.km1 = fitpar(8);
par.km2 = fitpar(9);

t1 = par.t;
phitest1 = par.phi;

Pcmb = zeros(length(t1),length(phitest1));
Pcdat = zeros(length(t1),length(phitest1));

n = NaN(length(t1),length(phitest1));
data = struct('x',n(:,:,[1 1 1 1 1]),'A',NaN,'L',n,'V',n,'lc',n,'vc',n,'Pc',n,'Pcdat',n);

for i = 1:length(phitest1)
    phi1 = phitest1(i);
    
    A1 = 0.125 / par.L0;
    par.L = @(t) par.L1 + A1 * cos(2*pi * (t - phi1));
    par.V = @(t) -2*pi * A1 * sin(2*pi * (t - phi1));
    switch par.model
        case 'lc'
            X0 = [par.L1+A1   0   0   0   1];
        case 'ls'
            X0 = [0   0   0   0   1];
        case {'old','old2'}
            X0 = [0   0   0   0   1];
    end
    
    %options for ode
    odeopt = odeset('RelTol',1e-5); %, 'OutputFcn',@odeplot);
    sol1 = ode23(@(t,x) odefcn(t,x,par), [0 t1(end)+1], X0, odeopt);
    
    x1 = deval(sol1,t1+1);
    
    [~,lc1,vc1,Pc1] = odefcn(t1',x1,par);
    x1 = x1';
    lc1 = lc1';
    vc1 = vc1';
    Pc1 = Pc1';
    
    Pcmb(:,i) = Pc1;
    Pcdat(:,i) = interp1(par.tdat, par.Pdat(:,i), t1) / par.P0;
    
    data.x(:,i,:) = permute(x1,[1 3 2]);
    data.A = A1;
    data.L(:,i) = par.L(t1);
    data.V(:,i) = par.V(t1);
    data.lc(:,i) = lc1;
    data.vc(:,i) = vc1;
    data.Pc(:,i) = Pc1;
    data.Pcdat(:,i) = Pcdat(:,i);
end
%dx = max(Pcmb) - max(Pcdat);
dx = flatten(Pcmb - Pcdat);
dx(isnan(dx)) = 0;
Pcmb = Pcmb(:);

fprintf('m = %g; b = %g; lc0 = %g; k1 = %g; k2 = %g; k30 = %g; k40 = %g; km1 = %g; km2 = %g    ->> sum(dx^2) = %f\n', ...
par.mm, par.b, par.lc0, par.k1, par.k2, ...
par.k30, par.k40, par.km1, par.km2, sum(dx.^2));

