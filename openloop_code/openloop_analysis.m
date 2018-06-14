%==========================================================================
%   This programm calculates the inherent motion
%   characteristics of the reduced F-16 model
%   - periodic: short period, phugoid, Dutch roll
%   - aperiodic: aperiodic roll, spiral
%
%   Author: Lukas
% 
%==========================================================================
clear; 

% import matrices
load redu_ss_openloop

% define trim point
v0 = 600; %[ft/s]

% set plotting/display parameters
pzmaps          = false;
characteristics = false;
simulations     = false;

% ------------------SS matrices--------------------------------------------

% LONGITUDINAL
states_lon    = {'theta' 'v' 'alpha' 'q'};
inputs_lon    = {'elevator'};
outputs_lon   = {'theta' 'v' 'alpha' 'q'};
SS_lon = ss(A_lon, B_lon, C_lon, D_lon,...
            'statename',states_lon,...
            'inputname',inputs_lon,...
            'outputname',outputs_lon);
        
% LATERAL
states_lat    = {'phi' 'beta' 'roll rate' 'yaw rate'};
inputs_lat    = {'aileron' 'rudder'};
outputs_lat   = {'phi' 'beta' 'roll rate' 'yaw rate'};
SS_lat = ss(A_lat, B_lat, C_lat, D_lat,...
            'statename',states_lat,...
            'inputname',inputs_lat,...
            'outputname',outputs_lat);

% -------------------pole-zero maps----------------------------------------
if pzmaps
    
    % LONGITUDINAL
    figure(1); 
    pzmap(SS_lon);
    print -depsc2 -r1200 Figures/pzmap_sym
    
    % LATERAL
    figure(2); 
    pzmap(SS_lat);
    print -depsc2 -r1200 Figures/pzmap_asym
    
end

% ------------------characteristics----------------------------------------
if characteristics
    
    % calculate natural frequencies and damping ratios
    [Wn_lon, zeta_lon, Po_lon] = damp(SS_lon);
    [Wn_lat, zeta_lat, Po_lat] = damp(SS_lat);
    
    % LONGITUDINAL
    
    % calculate period (PDF P.5 check!)
    P_lon = 2*pi()./(Wn_lon.*sqrt(1-zeta_lon.^2));
    % calculate time to half amplitude
    T_halfa_lon = log(2)./(Wn_lon.*zeta_lon);
    % combine characteristics into matrix
    chara_lon = [real(Po_lon),imag(Po_lon),Wn_lon,zeta_lon, P_lon, T_halfa_lon];
    
    % display matrix
    disp('Longitudinal:')
    disp('   Pole (real,imag)       Wn       zeta      P      T_halfa')
    disp('   --------------------------------------------------------')
    disp(chara_lon)
    
    % wrtie matrix to excel
    xlswrite('Tables/chara_lon',chara_lon);
    
    % LATERAL
    
    % calculate time constant
    tau_lat = 1./(zeta_lat.*Wn_lat);
    % calculate time to half amplitude
    T_halfa_lat = log(2)./(Wn_lat.*zeta_lat);
    % combine characteristics into matrix
    chara_lat = [real(Po_lat), imag(Po_lat), Wn_lat,zeta_lat, tau_lat, T_halfa_lat];
    
    % display matrix
    disp('Lateral:')
    disp('   Pole (real,imag)       Wn       zeta      tau      T_halfa')
    disp('   ----------------------------------------------------------')
    disp(chara_lat)
    
    % wrtie matrix to excel
    xlswrite('Tables/chara_lat',chara_lat);
    
end

% ---------------------simulation------------------------------------------

if simulations
    
    % figure parameters for printing
    simfig_width    = 1000;
    simfig_height   = 400;
    
    % lsim parameters for pulse-shaped inputs
    T = 10; % observaton period [s]
    dt = 0.01; % sampling time [s]
    
    % pulse parameters
    d_ail = 1; % aileron pulse duration [s]
    m_ail = rad2deg(0.025); % aileron pulse magnitude [deg]
    
    d_rud = 1; % rudder pulse duration [s]
    m_rud = rad2deg(0.025); % rudder pulse magnitude [deg]
    
    % define time axis
    t  = 0:dt:T; N = length(t);
    
    % create signals
    N_ail = d_ail/dt;
    N_rud = d_rud/dt;
    
    sig_ail  = [m_ail*ones(1,N_ail),zeros(1,N-N_ail)];
    sig_rud  = [m_rud*ones(1,N_rud),zeros(1,N-N_rud)];
    sig_null = zeros(1,N);
    
    
    % SHORT PERIOD
    [y_sp, t_sp] = step(-0.05*SS_lon,10); % elevator input: -0.05 [rad]
    theta_sp = y_sp(:,1);
    q_sp = y_sp(:,4);
    
    figure('pos',[100 100 simfig_width simfig_height/2])
    figure(1)
    plot(t_sp,q_sp)
    xlabel('t [sec]')
    ylabel('q [deg/s]')
    grid on
    print -depsc2 -r1200 Figures/short_period
    
    % PHUGOID
    [y_ph, t_ph] = step(-0.05*SS_lon,300); % elevator input: -0.05 [rad]
    v_ph = y_ph(:,2) + v0; % trim velocity added
    alpha_ph = y_ph(:,3);
    
    figure('pos',[100 100 simfig_width simfig_height])
    figure(2)
    subplot(2,1,1)
    plot(t_ph,v_ph)
    xlabel('t [sec]')
    ylabel('V [ft/sec]')
    grid on
    
    subplot(2,1,2)
    plot(t_ph,alpha_ph)
    xlabel('t [sec]')
    ylabel('\alpha [deg]')
    grid on
    print -depsc2 -r1200 Figures/phugoid

    
    % DUTCH ROLL
    % use rudder pulse
    sig_dr = [sig_null',sig_rud'];
    y_dr = lsim(SS_lat,sig_dr,t);
    p_dr = y_dr(:,3);
    r_dr = y_dr(:,4);
    
    figure('pos',[100 100 simfig_width simfig_height])
    figure(3)
    subplot(2,1,1)
    plot(t,p_dr)
    xlabel('t [sec]')
    ylabel('p [deg/sec]')
    grid on
    
    subplot(2,1,2)
    plot(t,r_dr)
    xlabel('t [sec]')
    ylabel('r [deg/s]')
    grid on
    print -depsc2 -r1200 Figures/dutch_roll
    
    % APERIODIC ROLL
    % use aileron pulse
    sig_ap = [sig_ail',sig_null'];
    y_ap = lsim(SS_lat,sig_ap,t);
    p_ap = y_ap(:,3);
    
    % (Plot combined with spiral)
    
    % SPIRAL
    % use aileron pulse
    
    sig_spir = [sig_ail',sig_null'];
    y_spir = lsim(SS_lat,sig_spir,t);
    phi_spir = y_spir(:,1);
    
    figure('pos',[100 100 simfig_width simfig_height])
    figure(4)
    subplot(2,1,1)
    plot(t,p_ap)
    xlabel('t [sec]')
    ylabel('p [deg/sec]')
    grid on
    
    subplot(2,1,2)
    plot(t,phi_spir)
    xlabel('t [sec]')
    ylabel('\phi [deg]')
    grid on
    print -depsc2 -r1200 Figures/aperoll_spiral
    
end




