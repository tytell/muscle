function d = check_jacobian(t,X,dx, statefcn,jacobianfcn)

Jn = zeros(length(X), length(X));

for i = 1:length(X)
    X1 = X;
    X1(i) = X1(i) + dx(i);
    Yp = statefcn(t,X1);
    
    X1 = X;
    X1(i) = X1(i) - dx(i);
    Ym = statefcn(t,X1);
    
    dY1 = (Yp - Ym) / (2*dx(i));
    
    Jn(:,i) = dY1;
end

Ja = jacobianfcn(t,X);

fprintf('Norm of difference in analytical and numerical Jacobian: %f\n', norm(Ja-Jn));
fprintf('Difference matrix:\n');
disp(Ja-Jn);
d = norm(Ja-Jn);

