% ------------------Sizing --------------------% 
clear; clc; close all; 

% ----------------- Solar Arrays ---------------------% 
% Input Data: 

mu = astroConstants(13); % [km^3 / s^2]
Re = astroConstants(23); %km Radius of Earth
R_orb = 7714.4278; % [km] AVISO
h = R_orb - Re; % [km] Altitude
n = sqrt(mu/ R_orb^3); % mean angular velocity (rad/s)
T_orb = (2*pi)/n; % Period (s)

 alpha_eclipse = rad2deg(acos(Re/(Re+h)));
% 
t_sunlit_perc = (180 + (2*alpha_eclipse))/(360);
Td = T_orb*t_sunlit_perc; % [s] spacecraft time in daylight 
t_eclipse_perc =  (180 - (2*alpha_eclipse))/(360);
Te= T_orb*t_eclipse_perc; % [s] spacecraft time in eclipse  

% Mission Requirements 
%Td= 77*60; % [s] spacecraft time in daylight 
%Te= 35*60; % [s] spacecraft time in eclipse 
Xd=0.85; % [-] line efficiency in daylight 
Xe=0.65; % [-] line efficiency in eclipse
Pd=540 + 0.2*540; % [W] power request in daylight + margin of 20% has to be added 
Pe=250 + 0.2*250; % [W] power request in eclipse + margin of 20% has to be added 
I0=1358; % [W/m^2] Sun irradiance 
T_life=5; % [y] mission duration in YEARS 
theta=deg2rad(23.44); % [rad] inclination angle between array surface normal and the Sun direction 

% Characteristics of solar arrays 
eps_BOL=0.094; % [-] efficiency at beginning of life for highest T=95 °C
dpy = 0.0275 ; % [/year] degradation solar array 
%p = 57; % [W/m^2] specific power [COMPUTED]
I_D= 0.88; % [-] inherent degradation 
rho_SA=1/0.55; % [m^2/kg] not sure.... surface density
V_cell=2.6; % [V] single solar cell voltage (for the moment the same of ex session)
V_sys=28; % [V] overall system voltage 
A_cell=23.61*(10^-4); % [m^2] single solar cell surface (was in cm^2)
A_real = 9.5; % [m^2] real surface area of SA
h_solar_cell=sqrt(A_cell); % [m] height of one solar cell considering that a solar cell is squared
%-- Power request for SA for maximum power demand-- 
P_SA=(Pe*Te/(Xe*Td) + Pd/Xd); %[W]

% --- specific power output ----- 
P_in = eps_BOL* I0; % [W/m^2]

% ---- SA specific power at BOL---- 
P_BOL= P_in * I_D * cos(theta); % [W/m^2]

% --- Lifetime degradation ---
L_life= (1-dpy)^T_life; % [year]

% --- SA specific power at EOL ---
P_EOL=L_life*P_BOL; % [W/m^2]

% --- SA surface --- 
A_SA=P_SA/P_EOL; % [m^2] surface area of the solar arrays
% adding one string to take into account the failure
A_SA=A_SA+12*h_solar_cell;

% --- SA mass --- 
m_SA= A_SA/rho_SA; % [kg] mass of the solar arrays needed

% --- REFINED SIZING ---
% --- cells needed ----
N=ceil(A_SA/A_cell); % [-] 

% --- cells in series needed ----
N_series=ceil(V_sys/V_cell); % [-]

% --- actual voltage --- 
V_real = N_series*V_cell; % [V] 

% ---- Actual number of cells --- 
N_real = ceil(N/N_series)*N_series; % [-]

% --- Real surface area --- 
A_SA_real = N_real * A_cell; % [m^2]
% --- Real surface area --- 
m_SA_real = A_SA_real/rho_SA; % [m^2]

A_SA_real_mission=9.5; % [m^2] real area of solar panels 
m_SA_real_mission=0.55*9.5; % [kg] mass of solar panels in real mission

% Realtive error computation
err_A_rel = abs(A_SA_real-A_SA_real_mission)/A_SA_real_mission; % [-] relative error on real area
err_m_rel = abs(m_SA_real-m_SA_real_mission)/m_SA_real_mission; % [-] relative error on real mass

%% ------------------ Batteries -----------------------%

% Input Data: 
% Mission Requirements 
T_R=Te/3600; % [h] time window in which the battery must provide power  
P_R=540 + 0.2*540; % [W] power required for the most critical mode 

% Characteristics of batteries
Em=118; % [Wh/kg] specific energy 
E_v=230; % [Wh/dm^3] energy density 
eta=0.65; % [-] line efficiency 
DOD=0.4; % [-] Depth of Discharge 
N_batt = 1; % [-] number of batteries 
V_syst_batt=32.4; % [V]  system voltage 
V_cell_batt= 3.6; % [V] single battery cell voltage
mu=0.8; % [-] package efficiency 
C_cell= 27; % [Ah] single cell capacity 

% --- Required capacity --- 
C=T_R*P_R/(DOD*N_batt*eta); % [Wh]

% --- mass of batteries ---
m_batt=C/Em; % [kg] 

% --- volume of batteries ---
Vol_batt=C/E_v*1e-3; % [dm^3]

% --- Cells in series needed --- 
N_series_batt= ceil(V_syst_batt/V_cell_batt); % [-]

% ---- Actual voltage needed --- 
V_real_batt = N_series_batt*V_cell_batt; % [V]

% ---- Capacity of a single string --- 
C_string = mu*C_cell*V_real_batt; % [Wh]

% --- strings to put in parallel --- 
N_parallel_batt= ceil(C/C_string); % [-]

% --- actual battery system capacity ---
C_real = N_parallel_batt*C_string; % [Wh]

% -- 
D = 0.053; %[m^3]
h = 0.185; %[m^3]
Vol_real = pi*(D/2)^2*h*27;
m_real = 0.81*27; %[kg]
err_vol = abs(Vol_batt-Vol_real)/Vol_real;
err_mass = abs(m_batt-m_real)/m_real;