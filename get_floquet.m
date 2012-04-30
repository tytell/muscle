function data = get_floquet(data,jfcn, nmodes, varargin)

opt.neigenvalues = [];
opt = parsevarargin(opt,varargin, 4);

N2 = 2*nmodes + 1;
P = size(data.x,2);

tph = (1:N2)'/N2 * data.per;
x = deval(data.sol, tph)';

A0 = zeros(N2,P,P);
for m = 1:N2
    A0(m,:,:) = jfcn(tph(m),x(m,:)');
end

A = zeros(N2*P, N2*P);
for i = 1:P
    ki = (1:N2) + (i-1)*N2;
    for j = 1:P
        kj = (1:N2) + (j-1)*N2;
        
        A(ki,kj) = diag(A0(:,i,j));
    end
end

Ginv = make_hb_gamma(nmodes,N2,P);
Omega = make_hb_omega(nmodes,P);

Ahat = Ginv \ A * Ginv;

F = 2*pi/data.per * Omega - Ahat;

[uf,mu] = eig(F);

mu = diag(mu);
good = abs(imag(mu)) < 0.99*2*pi/data.per;

mu = -mu(good);
uf = uf(:,good);

[~,ord] = sort(real(mu), 'descend');
mu = mu(ord);
uf = uf(:,ord);

if (~isempty(opt.neigenvalues))
    mu = mu(1:opt.neigenvalues);
    uf = uf(:,1:opt.neigenvalues);
end
neig = length(mu);

uf = reshape(uf,[N2 P neig]);

x = zeros(length(data.t),P,neig);
for i = 1:P
    x(:,:,i) = fourier_sin_cos(data.t,uf(:,:,i),data.per);
end    

data.jfcn = jfcn;
data.fx = x;
data.fexp = mu;
data.fmode = uf;

function Ginv = make_hb_gamma(N,M, P)

ph = 2*pi*(1:M)' ./ M;
k = [1:N; 1:N];

k = k(:)';

ph2 = ph * k;

Ginv0 = ones(M, 2*N+1);
Ginv0(:,2:2:end) = cos(ph2(:,1:2:end));
Ginv0(:,3:2:end) = sin(ph2(:,2:2:end));

bd = cell(1,P);
[bd{:}] = deal(Ginv0);

Ginv = blkdiag(bd{:});



function Omega = make_hb_omega(N, P)

Omega0 = zeros(2*N+1);

for i = 1:N
    Omega0(2*i,   2*i+1) = i;
    Omega0(2*i+1, 2*i  ) = -i;
end

bd = cell(1,P);
[bd{:}] = deal(Omega0);

Omega = blkdiag(bd{:});

