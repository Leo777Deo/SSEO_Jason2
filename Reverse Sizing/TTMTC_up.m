clear all;
clc;

sat.D=0.19; % [m] diameter of the satellite antenna [not verified: typical value up to 0.24] 
    GS.D=18;
% Input  
P_input=20*10^3; % [W] input power [not verified: guessed]
R= 4*1e3; % [bps] data rate for the uplink [verified: MCU 9]
D_rx=sat.D;
f = 2.04049*1e9; % [Hz] frequency of the carrier [verified: NOAA overview of Jason-2]
mu_amp=0.15; % [-] amplification constant  [not verified: guessed] 
eta= 0.1; % [deg] pointing accuracy [not verified: typical value]
BER=10^-7; % [-] Bit error rate [not verified: taken from ex lesson]
Eb_2_No_min=5.2; % [dB] minimum error per bit to signal density for the BER selected [not verified: guessed]
alpha_enc=2; % [-] [convolutional]
alpha_mod=1; % [-] [not verified: BPSK modulation from http://www.astronautix.com/p/proteus.html]

L_atm=3; % [dB] atmospheric losses 
%L_rain=3; % [dB] rain losses 
B=30*1e6; %[Hz] receiver bandwidth [not verified: taken from ex lesson]

fprintf("The uplink sizing: \n")
% [No, Eb_2_No, SNR_carrier,L_space,P_mod_loss] = TTMTC_function(P_input,R,D_rx,f,mu_amp,Eb_2_No_min,alpha_enc,alpha_mod,eta,L_atm,L_rain,B)

k = 1.380648*10^(-23); % [m^2.kg.s^-2.K^-1] Botzmann constant 
    c= 3*1e8; % [m/s] speed of light
    
    % Constants of the antennas
    % Antennas' characteristics 
    sat.mu_ant=0.75; % [-] efficiency of the satellite antenna [not verified: taken from spiral antenna document]
    GS.mu_ant=0.55; % [-] efficiency of the GS antenna [not verified]
    sat.D=0.19; % [m] diameter of the satellite antenna [not verified: typical value up to 0.24] 
    GS.D=18; % [m] diameter of the GS antenna [not verified]  Wallops Island in Germany 
    L_cables= -3; % [dB] cable losses [not verified: typical value on worst case scenario]
    Ts=273; % [K] sensor temperature (characteristic of the antenna) [not verified: typical value taken]

    % Parameters of the orbit 
    r_max=1336*1e3; % [m] largest distance between S/C and GS [verified: eoportal]

    %------Computation Part-------
    % Frequency of carrier
    lambda=c/f;
    
    % Design: Transmitter power selection
    P_tx=10*log10(300);
    R_real=R*alpha_enc/alpha_mod;
    
    % Gain of antennas
    G_ant_sat=10*log10((pi*sat.D*sqrt(sat.mu_ant)/lambda)^2); % [dB] satellite antenna gain 
    G_ant_GS=10*log10((pi*GS.D*sqrt(GS.mu_ant)/lambda)^2); % [dB] GS antenna gain
    G_tx=G_ant_GS;
    G_rx=G_ant_sat;
    
    % Computation of the losses 
    theta_rx=65.3*lambda/D_rx;
    L_point=-12*(eta/theta_rx)^2;
    L_space=20*log10(lambda/(4*pi*r_max));
    L_atm_tot=-L_atm;
    
    % Power computations 
    EIRP=P_tx+G_tx+L_cables
    P_rx=EIRP+G_rx+L_space+L_atm_tot+L_point;
    No=10*log10(k*Ts); % [dB] System noise density 
    
    % Energy per bit over Noise density
    Eb_2_No=P_rx - No -10*log10(R_real)
    
    fprintf("Is the sizing of Eb_2_No correct ? \n")
    Eb_2_No>Eb_2_No_min+3 % boolean to know if the sizing is correct 
    
    beta_mod=78; % [deg] modulation index (depends on receiver) [not verified: value to be assumed (from ex lesson)]
    
    P_mod_loss=20*log10(cosd(beta_mod));
    P_carrier=P_rx+P_mod_loss;

    fprintf("Is the carrier power positive ? \n")
    P_carrier>0

    % Signal to Noise ratio
    SNR_carrier= P_carrier - No - 10*log10(B)
    SNR_min=10; % [dB] minimum SNR depending on receiver [not verified: taken from ex lesson]
    
    % Output of the sizing 
    fprintf("Is the receiver capable of tracking the signal? ")
    SNR_carrier>3+SNR_min

