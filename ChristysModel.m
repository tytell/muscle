% Script that will calculate the force produced by the muscle model for a
% variety of phases. I tried to write it very clearly and comment it
% extensively, if something is unclear, let me know.


clear all
close all
clc


set(0,'DefaultAxesLineStyleOrder','-|-.|--|:','DefaultLineLineWidth',3)
set(0,'DefaultTextFontSize',18)
set(0,'DefaultAxesFontSize',18)

mycol = 'gmkymcrgbkymcrgbk';



%% Define the parameters


t = 1.0;                                    % final time of the simulation

dt = 2.5e-05;                               % time step
L0 = 2.94;                                  % base length

alpha_m = 0.8;                              % slope of alpha for v_c<0
alpha_p = 2.9;                              % slope of alpha for v_c>=0
alphamax = 1.8 ;                            % upper limit of alpha
alphamin = 0.0;                             % lower limit of alpha
cc = 2;                                     % calcium sites
ss = 6;                                     % SR sites
dm = 0.2802;                                % damping term
km1 = 17.5804                               % WDD parameter 1
km2 = 6.0156;                               % WDD parameter 2
lc0 = 0.9678;                               % lc at maximum tenatic force
ls0 = 0.0;                                  % initial value of stretch element
mm = 0.0542;                                % mass term
ka = [6.7281, 23.2794, 51.3537, 19.3801];   % mass action rate constant
mu0 = 1;                                    % intercept of mu
mu1 = 23.0;                                 % slope of mu
lambda_2 = -20.0;                           % lambda parameter

P0 = 60.86;                                 % P0 in thelma's simulations

my_phase = [0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];               % phases we have experimental data for


for j = 1:length(my_phase)
    phase = my_phase(j);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%              Initialize the vectors                                  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    my_t = 0:dt:t;                              % discretize time
    
    % % L = ones(size(my_t));                       % set the length of the segment
    % % V = zeros(size(my_t));                      % set the velocity of the segment
    L= (0.918+0.0425*cos(2*pi*(my_t-phase)));
    V= (-2*pi*0.0425*sin(2*pi*(my_t-phase)));
    
    
    [sizem,sizen] = size(my_t);                 % get the dimension of the time vector
    firstm = floor(sizem/2);                    % I used this to split it up for the first 5 cycles...
    firstn = floor(sizen/2);                    % ...and the last 5
    
    alpha = ones(size(my_t));                   % intialize alpha
    p = zeros(size(my_t));                      % initialize p
    pc= zeros(size(my_t));                      % initialize pc
    
    caf = zeros(size(my_t));                    % initialize caf
    ca = zeros(size(my_t));                     % initialize ca
    m = ones(size(my_t));                      % initialize m
    lc = zeros(size(my_t));                     % initialize lc
    lambda = ones(size(my_t));                  % initialize lambda
    vc = zeros(size(my_t));                     % initialize vc
    vs = zeros(size(my_t));                     % initialize vs
    ls = zeros(size(my_t));                     % initialize ls
    
    newmu = zeros(size(my_t));                  % initialize mu
    ksim = zeros(size(my_t));                   % initialize the activation tracker
    
    newk3 = zeros(size(my_t));                  % for k3 varying due to WDD
    newk4 = zeros(size(my_t));                  % for k4 varying due to WDD
    
    dm_dt = zeros(size(my_t));                  % these are for the ODES
    dcaf_dt = zeros(size(my_t));
    dca_dt = zeros(size(my_t));
    dpc_dt = zeros(size(my_t));
    dls_dt = zeros(size(my_t));
    dvs_dt = zeros(size(my_t));
    
    %% Start the simulation
    for k = 1:length(my_t)
        
        
        %% Check the activation
        if (my_t(k)-floor(my_t(k)))<0.36
            
            ksim(k) = 1;
            
            
        else
            
            ksim(k) = 0;
            
        end
        
        
        %% Update k3 and k4 (when WDD is added)
        newk3(k) = ka(3)/sqrt(m(k));
        newk4(k) = ka(4)*sqrt(m(k));
        
        %% Calculate the update for Ca and Caf for the next time step
        dca_dt(k) = (newk4(k)*caf(k) - newk3(k)*ca(k))*(1-caf(k)) ...
            + (ksim(k))*ka(1)*(cc-ca(k)-caf(k)) ...
            + (1-ksim(k))*(ka(2))*(ca(k))*(cc-ss-ca(k)-caf(k));
        
        dcaf_dt(k) = -(newk4(k)*caf(k)-newk3(k)*ca(k))*(1-caf(k));
        
        %% Calculate mu for the current time step
        newmu(k) = (mu0 + mu1*caf(k));
        
        
        %% Calculate lc and vc for this time step
        lc(k) = (L(k)-ls(k));
        vc(k) = (V(k)-vs(k));
        
        
        %% Determine alpha and lambda for the current time step
        if vc(k) >= 0
            
            alpha(k) = 1+alpha_p*vc(k);
            
        else
            
            alpha(k) = 1+alpha_m*vc(k);
            
        end
        
        lambda(k) = 1 + lambda_2*(lc(k)-lc0).^2;
        
        %% Bound alpha and lambda
        alpha(k) = max(alphamin,min(alphamax,alpha(k)));
        lambda(k) = max(0,min(lambda(k),1));
        
        %% Calculate the current Pc
        pc(k) = alpha(k)*lambda(k)*caf(k);
        
        
        %% WDD
        
        if vc(k) >= 0
            
            dm_dt(k) = -km2*(m(k)-1);
            
        else
            
            dm_dt(k) = -km1*pc(k)*vc(k);
            
        end
        
        %% Calculate the vs update
        dvs_dt(k) = 1./mm*(pc(k) - dm*vs(k) - newmu(k)*(ls(k)-ls0));
        
        
        
        %% find the values of ls, vs, ca, caf, and p for the next time step.
        ls(k+1) = ls(k) + dt*(vs(k));
        vs(k+1) = vs(k) + dt*(dvs_dt(k));
        ca(k+1) = ca(k) + dt*(dca_dt(k));
        caf(k+1) = caf(k) + dt*(dcaf_dt(k));
        m(k+1) = m(k) + dt*(dm_dt(k));
        p(k+1) = pc(k);
        
        
        
        
        
    end
    
    figure(1)
    hold all
    plot((j-1)+my_t,pc,'m')
    % %     axis([0 10 0 2])
    hold all
    
end


