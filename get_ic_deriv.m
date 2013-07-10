function icp = get_ic_deriv(odefcn, jfcn, ics)

[~,dfdxp] = jfcn(0,ics', ones(size(ics))');
isimp = any(dfdxp - eye(length(ics)) ~= 0);

icp = zeros(size(ics));
eqn = odefcn(0,ics',icp');
icp(~isimp) = -eqn(~isimp);

[impic,impicp] = decic(@(t,x,xp) odefcnreduced(t,x,xp, odefcn,ics',icp',isimp), 0, ...
    ics(isimp)',ones(sum(isimp),1), icp(isimp)',[]);
%ic(isimp) = impic;
icp(isimp) = impicp';


function eqn = odefcnreduced(t,x,xp, odefcn,ic0,icp0,isimp)

ic0(isimp) = x;
icp0(isimp) = xp;

eqn1 = odefcn(t,ic0,icp0);
eqn = eqn1(isimp);

