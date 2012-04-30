function test_vdp

pars.mu = 3;
pars.B = 3;
pars.T = 5;

[fp,per,data] = get_limit_cycle(@(t,x) odefcn(t,x,pars), 0.02, pars.T, [0 1], ...
    'Display','iter-detailed', 'fixedperiod',true, 'initialcycles',2, 'TolX',1e-8, 'RelTol',1e-6);

plot(data.x(:,1),data.x(:,2));

data = get_floquet(data,@(t,x) jfcn(t,x,pars), 50);

test_floquet(data, [0 0.25 0.5 0.75],0.1, 'fixedperiod',true, 'RelTol',1e-6);

function dx = odefcn(t,x,pars)

dx = [pars.mu * (x(1) - 1/3 * x(1)^3 - x(2)); ...
    1/pars.mu * (x(1) - pars.B * sin(2*pi/pars.T * t))];
    

function J = jfcn(t,x, pars)

J = [ (1-x(1)^2) * pars.mu, -pars.mu; ...
      1/pars.mu,             0  ];
  
