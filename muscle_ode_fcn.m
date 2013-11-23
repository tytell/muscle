function [dx,lc,vc,Pcval] = muscle_ode_fcn(t,x, par)

switch par.model
    case 'lc'
        lc = x(1,:);
        vc = x(2,:);        
        Ca = x(3,:);
        Caf = x(4,:);
        m = x(5,:);
        
        actval = par.act(t);
        Lval = par.L(t);
        Vval = par.V(t);
        
        ls = Lval - lc;
        
        gact = actval; %g(actval);

        Pcval = Pc(lc,vc,Caf, par);
        
        dm = par.km1*Pcval.*h(-vc, par) - par.km2*(m-1).*g(vc+0.5, par);
        
        k3 = par.k30 ./ sqrt(m);
        k4 = par.k40 .* sqrt(m);
        
        dCaf = (k3 .* Ca - k4 .* Caf) .* (1 - Caf);
        dCa = (k4 .* Caf - k3 .* Ca) .* (1 - Caf) + ...
            gact .* par.k1 .* (par.C - Ca - Caf) + ...
            (1 - gact) .* par.k2 .* Ca .* (par.C - par.S - Ca - Caf);
        dlc = vc;
        
        muval = mu(Caf, par);
        
        dvc = 1/par.mm * (-Pcval - par.b*vc + muval.*ls);
        
        dx = [dlc; dvc; dCa; dCaf; dm];
        
    case 'ls'
        ls = x(1,:);
        vs = x(2,:);
        Ca = x(3,:);
        Caf = x(4,:);
        m = x(5,:);
        
        actval = par.act(t);
        Lval = par.L(t);
        Vval = par.V(t);
        
        lc = Lval - ls;
        vc = Vval - vs;
        
        gact = actval; %g(actval);

        Pcval = Pc(lc,vc,Caf, par);
        
        dm = par.km1*Pcval.*h(-vc, par) - par.km2*(m-1).*g(vc+0.5, par);
        
        k3 = par.k30 ./ sqrt(m);
        k4 = par.k40 .* sqrt(m);
        
        dCaf = (k3 .* Ca - k4 .* Caf) .* (1 - Caf);
        dCa = (k4 .* Caf - k3 .* Ca) .* (1 - Caf) + ...
            gact .* par.k1 .* (par.C - Ca - Caf) + ...
            (1 - gact) .* par.k2 .* Ca .* (par.C - par.S - Ca - Caf);
        dls = vs;
        
        muval = mu(Caf, par);
        
        dvs = 1/par.mm * (Pcval - par.b*vs - muval.*ls);
        
        dx = [dls; dvs; dCa; dCaf; dm];
        
    case 'old'
        
        P = x(1,:);
        Ca = x(3,:);
        Caf = x(4,:);
        m = x(5,:);
        
        actval = par.act(t);
        Lval = par.L(t);
        Vval = par.V(t);
        
        gact = actval; %g(actval);

        k3 = par.k30 ./ sqrt(m);
        k4 = par.k40 .* sqrt(m);
        
        dCaf = (k3 .* Ca - k4 .* Caf) .* (1 - Caf);
        dCa = (k4 .* Caf - k3 .* Ca) .* (1 - Caf) + ...
            gact .* par.k1 .* (par.C - Ca - Caf) + ...
            (1 - gact) .* par.k2 .* Ca .* (par.C - par.S - Ca - Caf);
        
        muval = mu(Caf, par);

        lc = Lval - P./muval;
        vc_sign = muval.* Vval - par.k5*(Caf.*lambda(lc, par) - P) + ...
            P.*par.mu1./muval.*dCaf;
        if (vc_sign < 0)
            alpha1 = par.alpham;
        else
            alpha1 = par.alphap;
        end
        
        vc = vc_sign./(muval + par.k5.*Caf.*lambda(lc, par).*alpha1);

        Pcval = Pc(lc,vc,Caf, par);
        dP = par.k5*(Pcval-P);
        
        %dP = (lambda(lc) * Caf * (1 + alpha1*Vval + alpha1*mu1*P * dCaf/muval^2) - P) / ...
        %    (1/k5 + lambda(lc) * alpha1 * Caf/muval);
                
        dm = par.km1*Pcval.*h(-vc,par) - par.km2*(m-1).*g(vc+0.5,par);
        
        dx = [dP; zeros(size(dP)); dCa; dCaf; dm];
        
    case 'old2'
        
        P = x(1,:);
        Ca = x(3,:);
        Caf = x(4,:);
        m = x(5,:);
        
        actval = par.act(t);
        Lval = par.L(t);
        Vval = par.V(t);
        
        gact = actval; %g(actval);

        muval = mu(Caf, par);

        lc = Lval - P./muval;
        lambdaval = lambda(lc,par);

        k1 = par.k1 .* gact;
        k2 = par.k2 .* (1-gact);
        k3 = par.k30 ./ sqrt(m);
        k4 = par.k40 .* sqrt(m);
        
        xboth = (1-Caf) .* (k3.*Ca - k4.*Caf);
        dCaf = xboth;
        dCa = k1 .* (par.C - Ca - Caf) + ...
            k2 .* Ca .* (par.C - par.S - Ca - Caf) - xboth;
        
        vc_sign = muval .* Vval - par.k5*(Caf.*lambdaval-P) + ...
            P.*par.mu1./muval.*xboth; %for var.mu

        a1 = par.alpham*(vc_sign<=0) + par.alphap*(vc_sign>0);
        
        vc = vc_sign./(muval + par.k5.*Caf.*lambdaval.*a1);
        
        dW = -vc.*P;
        
        dm = (dW > 0).*(par.km1*dW)+ (dW<=0).*(par.km2*(1-m));
        
        alphaval = min(par.alphamax,max(0,1+a1.*vc));
        Pcval=Caf.*lambdaval.*alphaval;
        
        dP = par.k5*(Pcval-P);

        dx = [dP; dW; dCa; dCaf; dm];
end



function [hx,dhx] = h(x, par)

if (x > 10*par.s)
    hx = x;
    dhx = 1;
elseif (x < -10*par.s)
    hx = 0;
    dhx = 0;
else
    exs = exp(x/par.s);
    hx = par.s * log(1 + exs);
    if (nargout == 2)
        dhx = exs ./ (1 + exs);
    end
end

function gx = g(x, par)

gx = 1./(1+exp(-2*(x-0.5)/par.s));

function [hl,dl] = lambda(lc, par)

l0 = 1 + par.lambda2 * (lc - par.lc0).^2;
hl = l0;
hl(hl < 0) = 0;
% if (nargout == 1)
%     hl = h(l0, par);
% else
%     [hl,dhl] = h(l0, par);
%     dl = 2.*par.lambda2.*(lc - 1) .* dhl;
% end

function [x,dx] = alpha(vc, par)

if (nargout == 1)
    x = zeros(size(vc));
    x(vc >= 0) = 1 + par.alphap * vc(vc >= 0);
    x(vc < 0) = 1 + par.alpham * vc(vc < 0);
    x(vc > par.alphamax) = par.alphamax;
    %x0 = 1 + alphap * h(vc) - alpham * h(-vc);
    %x = h(x0);
    %x = alphamax - h(alphamax - x);
else
    [hvcp,dhvcp] = h(vc,par);
    [hvcm,dhvcm] = h(-vc,par);
    x0 = 1 + par.alphap * hvcp - par.alpham * hvcm;
    [x,dhx] = h(x0,par);

    dx = (alpham .* dhvcm + alphap .* dhvcp) .* ...
        dhx;
end

function x = mu(Caf, par)

x = par.mu0 + par.mu1*Caf;

function p = Pc(lc, vc, Caf, par)

p = lambda(lc, par) .* alpha(vc, par) .* Caf;
