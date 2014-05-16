function test_two_muscles

%optimized parameters:
%2: m = 0.0542; b = 0.2802; lc0 = 0.9678; k1 = 6.7281; k2 = 23.2794; k30 = 51.3537; k40 = 19.3801; km1 = 17.5804; km2 = 6.0156    ->> sum(dx^2) = 6.056118

filename = '~etytel01/Matlab/muscle/test_two_muscles.mat';
quiet = true;
doplot = false;

par.L0 = 2.94;                  % mm
par.Lis = 2.7;                  % mm

par.lc0 = 1;                    % nondimensional
par.lambda2 = -20;              % nondim
par.L1 = par.lc0 - 0.05;        % nondim

par.C = 2;
par.S = 6;

par.P0 = 67;                    % mN/mm^2

par.k1 = 6.7281;                % 1/s
par.k2 = 23.2794;               % 1/s
par.k30 = 51.3537;              % 1/s
par.k40 = 19.3801;              % 1/s
par.km1 = 17.5804;              % 1/s
par.km2 = 6.0156;               % 1/s

par.mm = 0.0542;                % arbitrary
par.b = 0.2802;                 % arbitrary

par.alpham = 0.8;               
par.alphap = 2.9;
par.alphamax = 1.8;

par.mu0 = 1;
par.mu1 = 23;

par.s = 0.05;

par.A = 0.125/par.L0;
par.phi = 0.1;
par.T = 1;
par.duty = 0.36;

par.M = 0.5;
par.zeta = 0.5;
par.omegar = 2*pi*1.5;

inputyn('','clearsaved');

dt = 0.005;

omegarvals = 2*pi* ([0.3:0.05:1 1.2 1.5 2]);
showfreq = [1 7 15 18];

ls1ind = 1;
vs1ind = 2;
Ca1ind = 3;
Caf1ind = 4;
m1ind = 5;
ls2ind = 6;
vs2ind = 7;
Ca2ind = 8;
Caf2ind = 9;
m2ind = 10;
Lind = 11;
Vind = 12;

X0 = [0   0   0   0    1    ...
      0   0   0   0    1    ...
      0   0];
check_jacobian(0,X0', 0.02*ones(12,1), @(t,x) odefcn(t,x,par), @(t,x) jfcn(t,x,par));

if doplot
    tinit = [0 15*par.T];
    odeopt = odeset('RelTol',1e-6); %, 'OutputFcn', @odeplot);
    [t,x] = ode45(@(t,x) odefcn(t,x,par), tinit, X0, odeopt);

    lcL = x(:,Lind) + par.L1 - x(:,ls1ind);
    lcR = par.L1 - x(:,Lind) - x(:,ls2ind);
    vcL = x(:,Vind) - x(:,vs1ind);
    vcR = -x(:,Vind) - x(:,vs2ind);

    Pcl = Pc(lcL, vcL, x(:,Caf1ind), par);
    Pcr = Pc(lcR, vcR, x(:,Caf2ind), par);
    figureseries('Time series');
    clf;
    xx = [0 par.duty par.duty 0; 0.5 par.duty+0.5 par.duty+0.5 0.5]';
    yy = [0 0 1 1; 0 0 1 1]';

    xx = repmat(xx,[1 15]);
    xx(:,1:2:end) = xx(:,1:2:end) + repmat(0:14,[4 1]);
    xx(:,2:2:end) = xx(:,2:2:end) + repmat(0:14,[4 1]);
    yy = repmat(yy,[1 15]);

    hax(1) = subplot(2,1,1);
    hold on;
    fill(xx(:,1:2:end),yy(:,1:2:end), [0.5 0.5 1]);
    fill(xx(:,2:2:end),yy(:,2:2:end), [0.5 1 0.5]);
    plot(t,Pcl, t,Pcr,'LineWidth',2);
    hold off;
    ylabel('Left muscle');
    axis tight;

    hax(2) = subplot(2,1,2);
    addmplot(t,x(:,Lind),'k-', t,x(:,Vind),'k:', 'LineWidth',2);
    axis tight;

    linkaxes(hax, 'x');
    set(hax, 'XLim',[12 15]);
    set(hax(1), 'YLim',[0 2]);

    labellines(hax(1), {'left','right'}, 'location',[13.12 13.62]);
    labellines(hax(2), {'L','V'});
    print('-dpdf','test_two_muscles-0.pdf');
end

if (~getvar('-file',filename,'freqdata') || (~quiet && ~inputyn('Use existing frequency data?', 'default',true)))
    freqdata = struct([]);
    
    progress(0,length(omegarvals),'**** Omegar test');
    for i = 1:length(omegarvals)
        par.omegar = omegarvals(i);
        
        X0 = [0   0   0   0    1    ...
              0   0   0   0    1    ...
              0   0];
        [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), 0.005, par.T, X0, ...
                'Display','final', 'fixedperiod',true, 'initialcycles',5, 'TolX',1e-8, 'RelTol',1e-6);
        
        data1.lc = par.L1 + data1.x(:,Lind)*[1 -1] - data1.x(:,[ls1ind ls2ind]);
        data1.vc = data1.x(:,Vind)*[1 -1] - data1.x(:,[vs1ind vs2ind]);

        data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,[Caf1ind Caf2ind]), par);

        data1.omegar = par.omegar;
        data1.zeta = par.zeta;
        data1 = get_floquet(data1,@(t,x) jfcn(t,x,par), 150);
        freqdata = makestructarray(freqdata,data1);
        
        progress(i);
    end
    
    putvar('-file',filename,'freqdata');
end

if doplot
    figureseries('Time series res freq');
    clf;
    xx = [0 par.duty par.duty 0; 0.5 par.duty+0.5 par.duty+0.5 0.5]';
    yy = [0 0 30 30; 0 0 30 30]';

    xx = repmat(xx,[1 2]);
    xx(:,1:2:end) = xx(:,1:2:end) + repmat(0:1,[4 1]);
    xx(:,2:2:end) = xx(:,2:2:end) + repmat(0:1,[4 1]);
    yy = repmat(yy,[1 2]);

    ex = [5 15 17];
    t = freqdata(ex(1)).t;
    Pc1 = cat(2,freqdata(ex).Pc);
    Pc1 = Pc1(:,1:2:end);
    Pc2 = freqdata(ex(2)).Pc(:,2);

    L1 = zeros(size(Pc1));
    for i = 1:3
        L1(:,i) = freqdata(ex(i)).x(:,Lind);
    end

    hax(1) = subplot(2,1,1);
    hold on;
    fill(xx(:,1:2:end),yy(:,1:2:end), [0.7 0.7 0.7], 'EdgeColor','none');
    fill(xx(:,2:2:end),yy(:,2:2:end), 'w', 'EdgeColor','k');
    hln = plot([t; t+1],[Pc1; Pc1],'LineWidth',2);
    hln(end+1) = addplot([t; t+1],[Pc2; Pc2], 'g--');
    hold off;
    axis tight;

    hax(2) = subplot(2,1,2);
    addplot([t; t+1],[L1; L1], 'LineWidth',2);
    axis tight;

    linkaxes(hax, 'x');
    labellines(hln,{'2','1','f_{act}/f_{res} = 0.67','right'}, 'location',[0.08 0.25 0.35, 0.8], ...
               'rotation',0)
    print('-dpdf','test_two_muscles-1.pdf');

    figureseries('Effect of resonant frequency');
    xx = cat(3, freqdata.x);
    Pcall = cat(3, freqdata.Pc);

    subplot(1,2,1);
    imagesc(freqdata(1).t, 1./(omegarvals/(2*pi)), squeeze(xx(:,9,:))');
    ylabel('Frequency ratio f_{act}/f_{res} (Hz)');
    xlabel('Time (sec)');
    hcol = colorbar;
    ylabel(hcol, 'L');

    subplot(1,2,2);
    imagesc(freqdata(1).t, 1./(omegarvals/(2*pi)), squeeze(Pcall(:,1,:))' / P0);
    ylabel('Frequency ratio f_{act}/f_{res} (Hz)');
    xlabel('Time (sec)');
    hcol = colorbar;
    ylabel(hcol, 'Force ratio P_{c,1} / P_0');
    print('-dpdf','test_two_muscles-2.pdf');

    fxx = cat(4, freqdata.fx);
    sgn1 = sign(fxx(1,1,1,:));
    fxx(:,:,1,:) = bsxfun(@times, fxx(:,:,1,:), sgn1);

    fexp = cat(2,freqdata.fexp);

    figureseries('Floquet exponents vs resonant frequency');
    clf;
    plot(1./(omegarvals/(2*pi)), log(0.5) ./ real(fexp(1:2,:))');
    xlabel('Frequency ratio f_{act}/f_{res} (Hz)');
    ylabel('t_{1/2} (sec)');
    title('Mode one time constants');
    print('-dpdf','test_two_muscles-3.pdf');

    figureseries('Floquet modes vs resonant frequency');
    clf;
    hln = -1*ones(10,1);
    for i = 1:4,
        subplot(2,2,i);
        j = showfreq(i);
        
        if (isreal(fxx(:,:,1,j)))
            hln(1:4) = plot(freqdata(j).t, fxx(:,1:4,1,j));
            hln(5:8) = addplot(freqdata(j).t, fxx(:,5:8,1,j),'--');
            hln(9:10) = addmplot(freqdata(j).t, fxx(:,9:10,1,j),'k-|k--','LineWidth',2);
        else
            plot(freqdata(j).t, real(fxx(:,:,1,j)));
            addplot(freqdata(j).t, imag(fxx(:,:,1,j)),'--');
        end
        xlabel('Time (s)');
        title(sprintf('f_{res} = %g',omegarvals(j)/(2*pi)));
    end
    legend(hln([1:4 9 10]),'lc','vc','Ca','Caf','L','V','Location','best');
    print('-dpdf','test_two_muscles-4.pdf');

    odeopt = odeset('RelTol',1e-6); %, 'OutputFcn', @odeoutput);
    par.omegar = 2*pi*0.5;
    i = find(omegarvals == par.omegar);
    ph = 0.2;
    j = first(freqdata(i).t >= ph);

    X0 = freqdata(i).x(j,:);
    X0 = X0 + fxx(j,:,1,i)*0.3;
    [t,x] = ode45(@(t,x) odefcn(t,x,par), freqdata(i).t(j) + [0 2], X0, odeopt);

    Pc1 = freqdata(i).Pc;
    Pcpert(:,1) = Pc(x(:,1), x(:,2), x(:,4));
    Pcpert(:,2) = Pc(x(:,5), x(:,6), x(:,8));

    figureseries('Position perturbation');
    hax = zeros(2,1);
    hax(1) = subplot(2,1,1);
    plot([freqdata(i).t; freqdata(i).t+1], [freqdata(i).x(:,9); freqdata(i).x(:,9)]-L0, 'k-', ...
         'LineWidth',1.5);
    addplot(t,x(:,9)-L0, 'LineWidth',2);

    ylabel('L (cm)');
    xtick labeloff;

    hax(2) = subplot(2,1,2);
    plot([freqdata(i).t; freqdata(i).t+1], [freqdata(i).Pc; freqdata(i).Pc], 'k-', ...
         'LineWidth',1.5);
    addmplot(t,Pcpert,'r-|r--', 'LineWidth',2);

    ylabel('P_c (mN)');
    xlabel('Time (sec)');
    %** Print figure
    print('-dpdf','test_two_muscles-5.pdf');
end

zetaold = par.zeta;
omegarold = omegarvals;

zetavals = [0.2 1 2 4];
omegarvals = 2*pi* ([0.3 0.5 0.8 1 1.2 1.5 2]);

if (~getvar('-file',filename,'dampdata') || (~quiet && ~inputyn('Use existing damping data?', 'default',true)))
    dampdata = struct([]);
    
    X0 = [0   0   0   0    1    ...
          0   0   0   0    1    ...
          0   0];
    
    a = 1;
    n = length(omegarvals) * length(zetavals);
    progress(0,n, '**** Daming test');
    for j = 1:length(zetavals)
        par.zeta = zetavals(j);
        for i = 1:length(omegarvals)
            par.omegar = omegarvals(i);
            
            fprintf('Zeta = %f, OmegaR = %f\n', par.zeta, par.omegar);
            
            [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), 0.005, par.T, X0, ...
                'Display','final', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
            data1.lc = par.L1 + data1.x(:,Lind)*[1 -1] - data1.x(:,[ls1ind ls2ind]);
            data1.vc = data1.x(:,Vind)*[1 -1] - data1.x(:,[vs1ind vs2ind]);

            data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,[Caf1ind Caf2ind]), par);
            
            data1 = get_floquet(data1,@(t,x) jfcn(t,x, par), 150);
            data1.zeta = par.zeta;
            data1.omegar = par.omegar;
            dampdata = makestructarray(dampdata,data1);
            
            a = a+1;
            progress(a);
        end
    end
    putvar('-file',filename,'dampdata');    
end

if doplot
    isomega = ismember(omegarold,omegarvals);
    islowzeta = zetavals < zetaold;

    dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)]);
    dampdata = makestructarray(dampdata(:,islowzeta), freqdata(isomega), dampdata(:,~islowzeta));
    dampdata = reshape(dampdata, [length(omegarvals) length(zetavals)+1]);
    fexp = catuneven(2,dampdata.fexp);
    fexp = reshape(fexp,[size(fexp,1) length(omegarvals) length(zetavals)+1]);

    for i = 1:2
        switch i
          case 1
            figureseries('Damping and resonance, mode 1');
          case 2
            figureseries('Damping and resonance, mode 2');
        end            
        clf;
        hax(i) = gca;
        plot(1./(omegarvals/(2*pi)), log(0.5) ./ real(squeeze(fexp(i,:,:))));
        xlabel('Frequency ratio f_{act}/f_{res} (Hz)');
        ylabel('t_{1/2} (sec)');
        title(sprintf('Mode %d time constants',i));
    end

    zetavalall = [zetavals(islowzeta) zetaold zetavals(~islowzeta)];
    lab = cell(size(zetavalall));
    for i = 1:length(zetavalall),
        lab{i} = sprintf('%g',zetavalall(i));
    end
    labellines(hax(1), lab,'location',-0.2);
    %** Print figure
    print('-dpdf','test_two_muscles-6.pdf');


    figureseries('Damping and floquet modes');
    clf;

    showdamp = [1 1; length(omegarvals) 1; 1 length(zetavalall); length(omegarvals) length(zetavalall)];
    for k = 1:4
        subplot(2,2,k);
        i = showdamp(k,1);
        j = showdamp(k,2);
        
        hln = zeros(6,1);
        hln(1:4) = plot(dampdata(i,j).t, dampdata(i,j).fx(:,1:4,1));
        addplot(dampdata(i,j).t, dampdata(i,j).fx(:,5:8,1),'--');
        hln(5:6) = addmplot(dampdata(i,j).t, dampdata(i,j).fx(:,9:10,1), 'k-|k--', 'LineWidth',2);
        xlabel('Time (sec)');
        
        title(sprintf('f_{act}/f_{res} = %g, \\zeta = %g, t_{1/2} = %g', ...
                      1./(omegarvals(i)/(2*pi)), zetavalall(j), log(0.5) ./ real(fexp(1,i,j))));
    end
    labellines(hln, {'lc','vc','Ca','Caf','L','V'});
    %** Print figure
    print('-dpdf','test_two_muscles-7.pdf');

    figureseries('Effect of the mode');
    clf;

    i = 1;
    j = length(zetavalall);
    plot(dampdata(i,j).t, dampdata(i,j).fx(:,1:4,1));
    addplot(dampdata(i,j).t, dampdata(i,j).fx(:,5:8,1),'--');
    addmplot(dampdata(i,j).t, dampdata(i,j).fx(:,9:10,1), 'k-|k--', 'LineWidth',2);
    xlabel('Time (sec)');
    ylabel('Deviation');

    labellines(hln, {'lc1','vc1','Ca1','Caf1','lc2','vc2','Ca2','Caf2','L','V'});
    %** Print figure
    print('-dpdf','test_two_muscles-8.pdf');
end

zetaoldvals = zetavals;
omegaroldvals = omegarvals;
dutycycleold = 0.36;

dutycyclevals = [0.2 0.5 0.7];
zetavals = [0.2 4];
omegarvals = 2*pi* ([0.3 2]);
if (~getvar('-file',filename,'dutydata') || (~quiet && ~inputyn('Use existing damping data?', 'default',true)))
    dutydata = struct([]);
    
    X0 = [0   0   0   0    1    ...
          0   0   0   0    1    ...
          0   0];
    
    a = 1;
    n = length(omegarvals) * length(zetavals) * length(dutycyclevals);
    progress(0,n, '**** Duty cycle tests');
    for k = 1:length(dutycyclevals)
        par.duty = dutycyclevals(k);
        
        for j = 1:length(zetavals)
            par.zeta = zetavals(j);
            for i = 1:length(omegarvals)
                par.omegar = omegarvals(i);
                fprintf('Duty = %g, Zeta = %g, OmegaR = %g\n', par.duty, par.zeta, par.omegar);

                [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), 0.005, par.T, X0, ...
                    'Display','final', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
                data1.lc = par.L1 + data1.x(:,Lind)*[1 -1] - data1.x(:,[ls1ind ls2ind]);
                data1.vc = data1.x(:,Vind)*[1 -1] - data1.x(:,[vs1ind vs2ind]);

                data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,[Caf1ind Caf2ind]), par);

                data1 = get_floquet(data1,@(t,x) jfcn(t,x, par), 150);
                
                data1.dutycycle = par.duty;
                data1.zeta = par.zeta;
                data1.omegar = par.omegar;
                dutydata = makestructarray(dutydata,data1);

                a = a+1;              
                progress(a);
            end
        end
    end
    putvar('-file',filename, 'dutydata');
end

if doplot
    isomega = ismember(omegaroldvals,omegarvals);
    iszeta = ismember(zetaoldvals,zetavals);
    islowduty = dutycyclevals < dutycycleold;

    dutydata = reshape(dutydata, [length(omegarvals) length(zetavals) length(dutycyclevals)]);
    dutydata = makestructarray(dutydata(:,:,islowduty), dampdata(isomega,iszeta), dutydata(:,:,~islowduty));
    dutydata = reshape(dutydata, [length(omegarvals) length(zetavals) length(dutycyclevals)+1]);

    dutyall = [dutycyclevals(islowduty) dutycycleold dutycyclevals(~islowduty)];

    figureseries('Force length vs duty cycle');
    xx = cat(3, dutydata.x);
    Lall = squeeze(xx(:,9,:));
    Lall = reshape(Lall,[size(Lall,1) 4 4]);

    Pcall = cat(3, dutydata.Pc);
    Pcall = squeeze(Pcall(:,1,:));
    Pcall = reshape(Pcall,[size(Pcall,1) 4 4]);

    for i = 1:4,
        hax(1) = subplot(2,2,i);
        mplot(Lall(:,:,i), Pcall(:,:,i), 'b-|g--|r-|c--','LineWidth',2);
    end
    %** Print figure
    print('-dpdf','test_two_muscles-9.pdf');

    
    fexp = catuneven(2,dutydata.fexp);
    fexp = reshape(fexp,[size(fexp,1) length(omegarvals) length(zetavals) length(dutyall)]);

    figureseries('Floquet exponents vs duty cycle');
    lab = {'0.2, 3', '0.2, 0.5', '4, 3', '4, 0.5'};
    clf;
    for i = 1:2
        hax(i) = subplot(2,1,i);
        plot(dutyall, log(0.5) ./ real(flatten(fexp(i,:,:,:),1:3)), 'o-');
        xlabel('Duty cycle');
        ylabel('t_{1/2} (sec)');
        title(sprintf('Mode %d time constants',i));
        legend(lab, 'Location','NE');
    end
    %** Print figure
    print('-dpdf','test_two_muscles-10.pdf');
end

vals = fullfact([2 2 2 1]);
islen = vals(:,1) == 2;
isvel = vals(:,2) == 2;
iswork = vals(:,3) == 2;
stiffval = vals(:,4);

par0 = par;
dutycyclevals = [0.2 0.36 0.5 0.7];
omegarvals = 2*pi* ([0.5 1 2]);
if (~getvar('-file',filename,'nonlindata') || (~quiet && ~inputyn('Use existing nonlinear data?', 'default',true)))
    nonlindata = struct([]);
    
    X0 = [0   0   0   0    1    ...
          0   0   0   0    1    ...
          0   0];
    par.zeta = zetaold;
    
    a = 1;
    n = length(dutycyclevals) * length(omegarvals) * length(islen);
    progress(0,n, '**** Nonlinear calculations');
    for k = 1:length(dutycyclevals)
        par.duty = dutycyclevals(k);
        for j = 1:length(omegarvals)
            par.omegar = omegarvals(j);
            for i = 1:length(islen)
                if (islen(i))
                    par.lambda2 = par0.lambda2;
                else
                    par.lambda2 = 0;
                end
                if (isvel(i))
                    par.alpham = par0.alpham;
                    par.alphap = par0.alphap;
                else
                    par.alpham = 0;
                    par.alphap = 0;
                end
                if (iswork(i))
                    par.km1 = par0.km1;
                else
                    par.km1 = 0;
                end
                switch stiffval(i)
                    case 1
                        par.mu0 = par0.mu0 + par0.mu1;
                        par.mu1 = 0;
                    case 2
                        par.mu0 = par0.mu0;
                        par.mu1 = 0;
                    case 3
                        par.mu0 = par0.mu0;
                        par.mu1 = par0.mu1;
                end
                
                fprintf('Duty = %g, OmegaR = %g, Nonlin = %d\n', par.duty, par.omegar, i);

                [~,~,data1] = get_limit_cycle(@(t,x) odefcn(t,x,par), 0.005, par.T, X0, ...
                    'Display','final', 'fixedperiod',true, 'initialcycles',10, 'TolX',1e-8, 'RelTol',1e-6);
                data1.lc = par.L1 + data1.x(:,Lind)*[1 -1] - data1.x(:,[ls1ind ls2ind]);
                data1.vc = data1.x(:,Vind)*[1 -1] - data1.x(:,[vs1ind vs2ind]);
                data1.Pc = Pc(data1.lc, data1.vc, data1.x(:,[Caf1ind Caf2ind]), par);

                data1 = get_floquet(data1,@(t,x) jfcn(t,x,par), 150);
                
                data1.dutycycle = par.duty;
                data1.zeta = par.zeta;
                data1.omegar = par.omegar;
                data1.lambda2 = par.lambda2;
                data1.islen = islen(i);
                data1.isvel = isvel(i);
                data1.iswork = iswork(i);
                data1.stiffval = stiffval(i);
                
                nonlindata = makestructarray(nonlindata,data1);
                a = a+1;
                
                progress(n)
            end
        end
    end
    putvar nonlindata;
end

if doplot
    nonlindata = reshape(nonlindata, [5 length(omegarvals) length(dutycyclevals)]);

    Pcall = cat(3, nonlindata.Pc);
    Pcall = reshape(Pcall,[size(Pcall,1) 2 4 3 4]);

    xx = cat(3, nonlindata.x);
    xx = reshape(xx, [size(xx,1) size(xx,2) size(nonlindata)]);
    tt = nonlindata(1).t;

    fexp = catuneven(2,nonlindata.fexp);
    fexp = reshape(fexp,[size(fexp,1) size(nonlindata)]);

    figureseries('Floquet exponents vs nonlinearity');
    lab = {'flat','FL','FV','FLV'};
    clf;
    for i = 1:2
        hax(i) = subplot(2,1,i);
        plot(dutycyclevals, log(0.5) ./ real(flatten(fexp(i,:,2,:),1:3)), 'o-');
        xlabel('Duty cycle');
        ylabel('t_{1/2} (sec)');
        title(sprintf('Mode %d time constants',i));
        legend(lab, 'Location','NE');
    end
    %** Print figure
    print('-dpdf','test_two_muscles-11.pdf');

    figureseries('Limit cycle vs nonlinearity');
    clf;
    subplot(2,1,1);
    plot(tt,squeeze(Pcall(:,1,:,2,2)));
    addplot(tt,squeeze(Pcall(:,2,:,2,2)),'--');

    subplot(2,1,2);
    plot(tt,squeeze(xx(:,9,:,2,2)));
    %** Print figure
    print('-dpdf','test_two_muscles-12.pdf');
    
    figureseries('Limit cycle vs duty cycle');
    clf;
    subplot(2,1,1);
    plot(tt,squeeze(Pcall(:,1,4,2,:)));

    subplot(2,1,2);
    plot(tt,squeeze(xx(:,9,4,2,:)));
    %** Print figure
    print('-dpdf','test_two_muscles-13.pdf');

    figureseries('Floquet exponents vs duty cycle2');
    clf;
    for i = 1:2
        hax(i) = subplot(2,1,i);
        plot(dutycyclevals, log(0.5) ./ real(flatten(fexp(i,4,:,:),1:3)), 'o-');
        xlabel('Duty cycle');
        ylabel('t_{1/2} (sec)');
        title(sprintf('Mode %d time constants',i));
    end
    %** Print figure
    print('-dpdf','test_two_muscles-14.pdf');
end



function actval = act(t, par)

t1 = mod(t,par.T);
actval(1,:) = double(t1 < par.duty);
actval(2,:) = double(((t1 >= 0.5) & (t1 < 0.5+par.duty)) | ((t1 >= 0) & (t1 < par.duty-0.5)));


function [hx,dhx] = h(x, par)

exs = exp(x/par.s);
hx = par.s * log(1 + exs);
if (nargout == 2)
    dhx = exs ./ (1 + exs);
end



function [gx,dgx] = g(x, par)

e2xs = exp(-2*x/par.s);
gx = 1./(1+e2xs);

if (nargout == 2)
    dgx = -2*e2xs ./ (1+e2xs);
end



function [hl,dl] = lambda(lc, par)

l0 = 1 + par.lambda2 * (lc - par.lc0).^2;
if (nargout == 1)
    hl = h(l0, par);
else
    [hl,dhl] = h(l0, par);
    dl = 2.*par.lambda2.*(lc - par.lc0) .* dhl;
end



function [x,dx] = alpha(vc, par)

if (nargout == 1)
    x0 = 1 + par.alphap * h(vc,par) - par.alpham * h(-vc,par);
    x = h(x0,par);
else
    [hvcp,dhvcp] = h(vc,par);
    [hvcm,dhvcm] = h(-vc,par);
    x0 = 1 + par.alphap * hvcp - par.alpham * hvcm;
    [x,dhx] = h(x0,par);

    dx = (par.alpham .* dhvcm + par.alphap .* dhvcp) .* ...
        dhx;
end

function p = Pc(lc, vc, Caf, par)

p = lambda(lc,par) .* alpha(vc,par) .* Caf;


function x = mu(Caf, par)

x = par.mu0 + par.mu1*Caf;


function dx = odefcn(t,x, par)

ls = x([1 6],:);
vs = x([2 7],:);
Ca = x([3 8],:);
Caf = x([4 9],:);
m = x([5 10],:);

L = x(11,:);
V = x(12,:);

actval = act(t,par);

gact = g(actval-0.5, par);

lc(1,:) = L + par.L1 - ls(1,:);
lc(2,:) = par.L1 - L - ls(2,:);
vc(1,:) = V - vs(1,:);
vc(2,:) = -V - vs(2,:);

Pcval = Pc(lc,vc,Caf, par);

dm = par.km1*Pcval .* h(-vc, par) - par.km2*(m-1).*g(vc,par);

k3 = par.k30 ./ sqrt(m);
k4 = par.k40 .* sqrt(m);

dCaf = (k3 .* Ca - k4 .* Caf) .* (1 - Caf);
dCa = (k4 .* Caf - k3 .* Ca) .* (1 - Caf) + ...
    gact .* par.k1 .* (par.C - Ca - Caf) + ...
    (1 - gact) .* par.k2 .* Ca .* (par.C - par.S - Ca - Caf);
dls = vs;

muval = mu(Caf, par);

dvs = 1/par.mm * (Pcval - par.b*vs - muval.*ls);

dL = V;
dV = 1/par.M * (-muval(1,:) .* ls(1,:) + par.b * vs(1,:) + ...
             muval(2,:) .* ls(2,:) - par.b * vs(2,:) ...
             - 2 * par.zeta * par.omegar * V - par.omegar^2 * L);

dx = [dls(1,:); dvs(1,:); dCa(1,:); dCaf(1,:); dm(1,:); ...
    dls(2,:); dvs(2,:); dCa(2,:); dCaf(2,:); dm(2,:); ...
    dL; dV];




function J = jfcn(t,x, par)

ls = reshape(x([1 6]), [1 1 2]);
vs = reshape(x([2 7]), [1 1 2]);
Ca = reshape(x([3 8]), [1 1 2]);
Caf = reshape(x([4 9]), [1 1 2]);
m = reshape(x([5 10]), [1 1 2]);
L = x(11);
V = x(12);

lc(1,1,1) = L + par.L1 - ls(1,1,1);
lc(1,1,2) = par.L1 - L - ls(1,1,2);
vc(1,1,1) = V - vs(1,1,1);
vc(1,1,2) = -V - vs(1,1,2);


[lambdaval,dlambdaval] = lambda(lc,par);
[alphaval,dalphaval] = alpha(vc,par);
actval = act(t, par);

gact = reshape(g(actval-0.5, par), [1 1 2]);

[hvc,dhvc] = h(vc, par);
[gvc,dgvc] = g(vc, par);

J = zeros(12,12);

%each muscle block has an independent 5x5 block in the Jacobian
%that describes its independent evolution
zero = zeros(1,1,2);
one = ones(1,1,2);
Jmusc = [ ...
    zero, one, zero, zero, zero; ...
    ...
    1/par.mm .* (-par.mu0 - Caf .* par.mu1 - Caf .* alphaval .* dlambdaval), ...
    1/par.mm .* (-par.b - Caf .* lambdaval .* dalphaval), ...
    zero, ...
    1/par.mm .* (-par.mu1 .* ls + lambdaval .* alphaval), ...
    zero; ...
    ...
    zero, ...
    zero, ...
    -(1 - Caf) .* par.k30 ./ sqrt(m) - Ca .* par.k2 .* (1 - gact) + ...
       par.k2 .*(par.C - par.S - Ca - Caf) .* (1 - gact) - par.k1 .* gact, ...
    Ca .* par.k30 ./ sqrt(m) + (1 - Caf).*par.k40.*sqrt(m) - Caf.*par.k40.*sqrt(m) - ...
       Ca.*par.k2 .* (1 - gact) - par.k1 .* gact, ...
    (1 - Caf) .* (Ca .* par.k30 ./ (2 * m.^1.5) + Caf .* par.k40 ./ (2*sqrt(m))); ...
    ...
    zero, ...
    zero, ...
    (1 - Caf) .* par.k30 ./ sqrt(m), ...
    -Ca .* par.k30 ./ sqrt(m) - (1 - Caf) .* par.k40 .* sqrt(m) + ...
       Caf .* par.k40 .* sqrt(m), ...
    (1 - Caf) .* (-Ca .* par.k30 ./ (2*m.^1.5) - Caf .* par.k40 ./ (2*sqrt(m))); ...
    ...
    -Caf .* par.km1 .* alphaval .* hvc .* dlambdaval, ...
    -Caf .* par.km1 .* hvc .* lambdaval .* dalphaval + ...
       par.km2 .* (m - 1) .* dgvc + Caf .* par.km1 .* alphaval .* lambdaval .* dhvc, ...
    zero, ...
    par.km1 .* alphaval .* hvc .* lambdaval, ...
    -par.km2 .* gvc
    ];

J(1:5,1:5) = Jmusc(:,:,1);
J(6:10,6:10) = Jmusc(:,:,2);


%spring mass block
J(11:12,11:12) = [0, 1; ...
    -par.omegar^2/par.M, -2*par.zeta*par.omegar/par.M];

%coupling between muscles and spring mass
J(12,1:10) = [1/par.M .* (-par.mu0 - Caf(1) .* par.mu1),   par.b/par.M, ...
    0, -par.mu1 .* ls(1) ./ par.M, 0, ...
    1/par.M .* (par.mu0 + Caf(2) .* par.mu1),   -par.b/par.M, ...
    0, -par.mu1 .* ls(2) ./ par.M, 0];

Jcouple = [zero, zero; ...
    ...
    Caf.*alphaval.*dlambdaval / par.mm,   ...
    Caf.*lambdaval.*dalphaval / par.mm; ...
    ...
    zero, zero; ...
    ...
    zero, zero; ...
    ...
    Caf .* par.km1 .* alphaval .* hvc .* dlambdaval,   ...
    Caf .* par.km1 .* hvc .* lambdaval .* dalphaval - ...
       par.km2 * (1 - m) .* dgvc - Caf .* par.km1 .* alphaval .* lambdaval .* dhvc(1)];

J(1:5,11:12) = Jcouple(:,:,1);
J(6:10,11:12) = Jcouple(:,:,2);

    



function stat = odeoutput(t,y, flag)

if (isempty(flag) && ~isempty(t) && (mod(t(1),1) < 0.0001))
    fprintf('t = %g\n', t(1));
end
stat = 0;

  

