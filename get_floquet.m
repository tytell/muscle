function data = get_floquet(data,jfcn, nmodes, varargin)

opt.odetype = 'direct';
opt = parsevarargin(opt, varargin, 4);

N2 = 2*nmodes + 1;
P = size(data.x,2);

tph = (1:N2)'/N2 * data.per;
[x,xp] = deval(data.sol, tph);
x = x';
xp = xp';

A0 = zeros(N2,P,P);
switch opt.odetype
    case 'direct'
        for m = 1:N2
            A0(m,:,:) = jfcn(tph(m),x(m,:)');
        end
    case 'implicit'
        for m = 1:N2
            [A1,B1] = jfcn(tph(m), x(m,:)', xp(m,:)');
            A0(m,:,:) = -B1 \ A1;
        end
end

A = sparse(N2*P, N2*P);
for i = 1:P
    ki = (1:N2) + (i-1)*N2;
    for j = 1:P
        kj = (1:N2) + (j-1)*N2;
        
        A = A + sparse(ki,kj, A0(:,i,j), N2*P,N2*P);
    end
end

Ginv = make_hb_gamma(nmodes,N2,P);
Omega = make_hb_omega(nmodes,P);

Ahat = Ginv \ A * Ginv;

F = Ahat - 2*pi/data.per * Omega;

%[uf,mu] = eig(F);
%[uf,mu] = eigs(F, P, 'si');
[uf,mu] = floquet_eigs(F, P, data.per);

uf = reshape(uf,[N2 P P]);

x = NaN(length(data.t),P,P);
for i = 1:P
    if ~isnan(mu(i))
        x(:,:,i) = fourier_sin_cos(data.t,uf(:,:,i),data.per);
    end
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
Ginv = sparse(Ginv);



function Omega = make_hb_omega(N, P)

a = 1:N;
i = [2*a; 2*a+1];
j = [2*a+1; 2*a];
s = [a; -a];

i = repmat(i,[1 1 P]);
j = repmat(j,[1 1 P]);
s = repmat(s,[1 1 P]);

i = bsxfun(@plus, i, reshape((2*N+1)*(0:P-1),[1 1 P]));
j = bsxfun(@plus, j, reshape((2*N+1)*(0:P-1),[1 1 P]));

Omega = sparse(i(:),j(:),s(:));

