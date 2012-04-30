function test_simple

[fp,per,data] = get_limit_cycle(@odefcn, 0.1, 10, [1.1 0], 'Display','iter-detailed');

plot(data.x(:,1),data.x(:,2));

data = get_floquet(data,@jfcn, 20);

function dx = odefcn(t,x)

r = sqrt(sum(x.^2));
dx = [-x(1) * (r-1) - x(2); -x(2) * (r-1) + x(1)];


function J = jfcn(t,x)

r = sqrt(sum(x.^2));
J = [ 1 - x(1)^2/r - r,      - (x(1)*x(2))/r - 1; ...
      1 - (x(1)*x(2))/r,     1 - x(2)^2/r - r];
  