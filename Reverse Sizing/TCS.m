%% Thermal Control System
clc;
clear all;
% Mono-nodal approach 
mu = astroConstants(13); % [km^3 / s^2]
m_sc_dry = 525; % kg from eoportal
m_TCS = 525*0.05; % kg mass of TCS
x_sc=1; %m  [presskit]
y_sc=1; %m
z_sc=3.7; %m
S_1 = x_sc*y_sc;
S_2 = x_sc*z_sc;
S_3 = S_2;
A_tot = (S_1*2) + (2*S_2) + (2*S_3); % [m^2]
alpha = 0.23;    % Electroplated Gold from NASA 
epsilon = 0.03;  % Electroplated Gold from NASA
sigma = 5.67*10^-8;  %Boltzmann constant
%T_Earth = 273.15 + 15; %[k] NASA Solar System Temperatures: Earth at 15 °C
T_Earth = 255; %[k] NASA Solar System Temperatures: Earth at 15 °C
T_DS = 3; %[k]
eps_rad = 0.76; %excercise session
Re = astroConstants(23); %km Radius of Earth
R_orb = 7714.4278; % [km] AVISO
h = R_orb - Re; % [km] Altitude
n = sqrt(mu/ R_orb^3); % mean angular velocity (rad/s)
T_orb = (2*pi)/n; % Period (s)

alpha_eclipse = rad2deg(acos(Re/(Re+h)));

t_sunlit_perc = (180 + (2*alpha_eclipse))/(360);
t_sunlit = T_orb*t_sunlit_perc;
t_eclipse_perc =  (180 - (2*alpha_eclipse))/(360);
t_eclipse = T_orb*t_eclipse_perc;

Tmax = 50; % [°C] max limit for electronics

Tmin = -20; % [°C]min limit for electronics


% Internal flux
Q_int_min = 180; % [w] No Payload + SHM, from MUP 3 and Jason series archive report 180
Q_int_max = 310 + 147; % [W] 300 for just platform and 147 for Payload only;Jason series archive report 211

% Sun flux
q_sun_sc = 1358;  %from attitude (at our altitude) used also in ADCS
Q_sun = S_1*alpha*q_sun_sc; % [W] S_1 as the Sun hits the spacecraft on the "head"

% Albedo
a = 0.39; % [-] Assumed from SSEO book between 0.31-0.39, 0.39 is worst case
q_alb = q_sun_sc*a; % [W/m^2]
Q_albedo = S_2*alpha*q_alb; % [W] S_2 as the Nadir-pointing face is the largest one and faces Earth

% Infrared 
eps_E = 0.95; % [-] emissivity of Earth
q_IR = eps_E * sigma * T_Earth^4; % [W/m^2]
Q_IR = S_2 * epsilon * q_IR; % [W] S_2 as the Nadir-pointing face is the largest one and faces Earth

% T spacecraft
A_DS_hot = S_1 + (3*S_2);% [m^2] The faces facing DS are the remaining ones, so bottom and 3 on the sides that are equal
T_sc_hot = ((Q_int_max + Q_sun + Q_albedo + Q_IR) / (sigma*epsilon*A_DS_hot) + T_DS^4)^(1/4); %[k]
T_sc_hot_cel = T_sc_hot -273.15 % [°C]

% Radiator area
T_max = -15+Tmax+273.15; % [K] assumed value of electronics (esercitatore ha detto così)

%A_rad = (Q_int_max + Q_sun + Q_albedo + Q_IR)/(sigma*eps_rad*(T_max^4-T_DS^4)) - epsilon*A_DS_hot/eps_rad; % [m^2] Minimum area to have T_sc < T_max
A_rad = (Q_int_max + Q_sun + Q_albedo + Q_IR - (sigma*epsilon*A_DS_hot*(T_max^4-T_DS^4)) )/ (sigma*(T_max^4-T_DS^4)*(eps_rad - epsilon)); % [m^2]
A_rad_perc = A_rad/A_tot;
T_sc_hot_check = ((Q_int_max + Q_sun + Q_albedo + Q_IR)/((sigma*epsilon*(A_DS_hot-A_rad))+(sigma*eps_rad*A_rad)))^(1/4);
T_sc_hot_check_cel = T_sc_hot_check -273.15 % [°C] Must be 50 °C :)

%% Cold case
A_DS_cold= (2*S_1) + (3*S_2);% [m^2] The faces facing DS are the remaining ones, so bottom and 3 on the sides that are equal
%T_sc_cold = ((Q_int_min + Q_IR) / (sigma*epsilon*A_DS_cold) + T_DS^4)^(1/4);  %[K]
T_sc_cold = ((Q_int_min + Q_IR) / ((sigma*epsilon*(A_DS_cold-A_rad))+(sigma*eps_rad*A_rad)) + T_DS^4)^(1/4);  %[K]
T_sc_cold_cel = T_sc_cold - 273.15 % [°C]
% Heaters 
T_min = Tmin+15+273.15; %[K] assumed value of electronics (esercitatore ha detto così)

%Q_heaters = sigma*epsilon*A_DS_cold*T_min^4 - Q_IR - Q_int_min; % [W] Heat of heaters
Q_heaters = ((sigma*epsilon*(A_DS_cold-A_rad))+(sigma*eps_rad*A_rad))*T_min^4 - Q_IR - Q_int_min; % [W] Heat of heaters
%T_sc_cold_check = ((Q_int_min + Q_IR + Q_heaters) / (sigma*epsilon*A_DS_cold) + T_DS^4)^(1/4);  %[K]
T_sc_cold_check = ((Q_int_min + Q_IR + Q_heaters) /((sigma*epsilon*(A_DS_cold-A_rad))+(sigma*eps_rad*A_rad)) + T_DS^4)^(1/4);  %[K]
T_sc_cold_check_cel = T_sc_cold_check -273.15 % [°C] Must be -20 °C and it is :)

%A_rad_cold = (Q_int_min + Q_IR)/(sigma*eps_rad*(T_min^4-T_DS^4)) - epsilon*(A_DS+A_rad)/eps_rad % [m^2] Minimum area to have T_sc < T_max

%% Panels: FLAT PLATE APPROX


b_p=1.490; %m MUP3
l_p=3.2885906; % m MUP3

eta_pan = 0.3;

A_panel = b_p*l_p; % m^2

eps_pan_hot = 0.82; % MUP 
alpha_pan_hot = 0.85; % MUP

eps_pan_cold = 0.7; % MUP
alpha_pan_cold = 0.92; % MUP

Q_int_pan = 580/2; % eta_pan*Q_sun_pan

Q_sun_pan = q_sun_sc*alpha_pan_hot*A_panel; %W

Q_alb_pan = q_alb*alpha_pan_hot*A_panel;  %W

%Q_IR_pan_hot = A_panel*eps_pan_hot*q_IR; %W

Q_IR_pan_cold = A_panel*eps_pan_cold*q_IR; %W

Q_sum_hot = Q_sun_pan + Q_alb_pan + Q_IR_pan_cold; % W
Q_sum_cold = Q_IR_pan_cold; % W

T_pan_hot = ((Q_sum_hot - (eta_pan*Q_sun_pan))/(sigma*(0.5*eps_pan_hot + 0.5*eps_pan_cold)*2*A_panel))^(1/4); % K
T_pan_hot_cel = T_pan_hot-273.15 % °C

T_pan_cold = ((Q_IR_pan_cold)/(sigma*(0.5*eps_pan_hot + 0.5*eps_pan_cold)*2*A_panel))^(1/4); % K
T_pan_cold_cel = T_pan_cold-273.15 % °C

 %% SPHERE s/c
% 
% A_sphere = A_tot; % [m^2] Area of equivalent sphere = Largest area of Jason-2
% r_sphere = sqrt(A_sphere/(4*pi)); % [m] Radius of equivalent sphere s/c
% A_cross = pi*r_sphere^2; % [m^2] Cross sectional area of equivalent sphere spacecraft
% 
% F = 0.5*(1 - ((sqrt(((h/Re)^2)+(2*(h/Re))))/(1+(h/Re)))); %[-] View Factor
% 
% % Internal flux
% Q_int_min = 300; % [w] No Payload, from MUP 3 and Jason series archive report
% Q_int_max = Q_int_min + 147; % [W] Jason series archive report
% 
% % Sun flux
% q_sun_sc = 1358;  %from attitude (at our altitude) used also in ADCS
% Q_sun = A_cross*alpha*q_sun_sc; % [W] S_1 as the Sun hits the spacecraft on the "head"
% 
% % Albedo
% a = 0.39; % [-] Assumed from SSEO book between 0.31-0.39, 0.39 is worst case
% q_alb = q_sun_sc*a*(Re/R_orb)^2; % [W/m^2]
% Q_albedo = A_sphere*alpha*q_alb*F; % [W] S_2 as the Nadir-pointing face is the largest one and faces Earth
% 
% % Infrared 
% eps_E = 0.95; % [-] emissivity of Earth
% q_IR = ((Re/R_orb)^2)*eps_E*sigma * T_Earth^4; % [W/m^2]
% Q_IR = A_sphere*F*epsilon*q_IR; % [W] S_2 as the Nadir-pointing face is the largest one and faces Earth
% 
% % T spacecraft
% A_DS_hot = A_sphere;
% T_sc_hot = ((Q_int_max + Q_sun + Q_albedo + Q_IR) / (sigma*epsilon*A_DS_hot) + T_DS^4)^(1/4); %[k]
% T_sc_hot_cel = T_sc_hot -273.15 % [°C]
% 
% % Radiator area
% T_max = 50+273.15; % [K] assumed value of electronics (esercitatore ha detto così)
% 
% %A_rad = (Q_int_max + Q_sun + Q_albedo + Q_IR)/(sigma*eps_rad*(T_max^4-T_DS^4)) - epsilon*A_DS_hot/eps_rad; % [m^2] Minimum area to have T_sc < T_max
% A_rad = (Q_int_max + Q_sun + Q_albedo + Q_IR - (sigma*epsilon*A_DS_hot*(T_max^4-T_DS^4)) )/ (sigma*(T_max^4-T_DS^4)*(eps_rad - epsilon)); % [m^2]
% 
% T_sc_hot_check = ((Q_int_max + Q_sun + Q_albedo + Q_IR)/((sigma*epsilon*(A_DS_hot-A_rad))+(sigma*eps_rad*A_rad)))^(1/4);
% T_sc_hot_check_cel = T_sc_hot_check -273.15 % [°C] Must be 50 °C :)
% 
% %% Cold case
% A_DS_cold= A_sphere;
% T_sc_cold = ((Q_int_min + Q_IR) / (sigma*epsilon*A_DS_cold) + T_DS^4)^(1/4);  %[K]
% T_sc_cold_cel = T_sc_cold - 273.15 % [°C]
% % Heaters 
% T_min = -20 + 273.15; %[K] assumed value of electronics (esercitatore ha detto così)
% 
% Q_heaters = sigma*epsilon*A_DS_cold*T_min^4 - Q_IR - Q_int_min; % [W] Heat of heaters
% 
% T_sc_cold_check = ((Q_int_min + Q_IR + Q_heaters) / (sigma*epsilon*A_DS_cold) + T_DS^4)^(1/4);  %[K]
% T_sc_cold_check_cel = T_sc_cold_check -273.15 % [°C] Must be -20 °C and it is :)
% 
% %A_rad_cold = (Q_int_min + Q_IR)/(sigma*eps_rad*(T_min^4-T_DS^4)) - epsilon*(A_DS+A_rad)/eps_rad % [m^2] Minimum area to have T_sc < T_max
