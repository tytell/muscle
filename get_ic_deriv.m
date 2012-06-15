function icp = get_ic_deriv(odefcn, jfcn, ics)

[~,dfdxp] = jfcn(0,ics', ones(size(ics))');
isimp = any(dfdxp - eye(length(ics)) ~= 0);

icp = zeros(size(ics));
eqn = odefcn(0,ics',icp');
icp(~isimp) = -eqn(~isimp);

[~,icp] = decic(odefcn, 0, ics',zeros(size(ics))', icp',~isimp);
icp = icp';

