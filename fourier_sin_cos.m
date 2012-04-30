function y = fourier_sin_cos(t, coefs,per)

N2 = size(coefs,1);

N = (N2-1)/2;

ph = 2*pi * t(:) / per;

k = [1:N; 1:N];

k = k(:)';

ph2 = ph * k;

Ginv = ones(size(ph2,1),2*N+1);
Ginv(:,2:2:end) = cos(ph2(:,1:2:end));
Ginv(:,3:2:end) = sin(ph2(:,2:2:end));

y = Ginv * coefs;

