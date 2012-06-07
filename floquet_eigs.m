function [uf,mu] = floquet_eigs(F,P, per, varargin)

opt.minfrac = 0.05;
opt.mult = 1.5;
opt = parsevarargin(opt, varargin, 4);

minexp = log(0.5) / (opt.minfrac * per);

guess = 0;
lasteig = 1;
n = 0;
mu = NaN(P,1);
uf = NaN(size(F,1),P);
while (n < P)
    [uf1, mu1] = eigs(F,1, guess);
    
    if ((abs(imag(mu1)) <= pi/per) && (abs(mu1 - lasteig)/abs(lasteig) > 0.001))
        if (~isreal(mu1))
            mu(n+(1:2)) = [mu1 conj(mu1)];
            uf(:,n+(1:2)) = [uf1 conj(uf1)];
            n = n+2;
        else
            mu(n+1) = mu1;
            uf(:,n+1) = uf1;
            n = n+1;
        end
        guess = real(mu1);
        lasteig = mu1;
    end
    
    if (real(mu1) < minexp)
        break;
    else
        guess = opt.mult * guess;
    end
end

        
        

