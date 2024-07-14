clear all;
clc;


Mass_sc=553; %kg
x_sc=1; %m  [presskit]
y_sc=1; %m
z_sc=3.7; %m

% x_sc=0.933; %m
% y_sc=0.839; %m
% z_sc=3.7; %m


Volume_sc=x_sc*y_sc*z_sc; % m^3

% Inerita matrix
Ix=(1/12) * Mass_sc *(y_sc^2 + z_sc^2); %kg*m^2
Iy=(1/12) * Mass_sc *(x_sc^2 + z_sc^2); %kg*m^2
Iz=(1/12) * Mass_sc *(y_sc^2 + x_sc^2); %kg*m^2
I=[Ix;Iy;Iz];
I=diag(I);
I_inv=I\eye(3,3);

% Faces  
r1=[y_sc/2;0;0]; % m
r2=[0;x_sc/2;0];
r3=[0;0;z_sc/2];
r4=[-y_sc/2;0;0];
r5=[0;-x_sc/2;0];
r6=[0;0;-z_sc/2];

A1=x_sc*z_sc; %m^2
A2=y_sc*z_sc;
A3=x_sc*y_sc;
A4=A1;
A5=A2;
A6=A3;

% Panels (m)  [proteus datasheet]

b_p=1.490; %m MUP3
l_p=3.2885906; % m MUP3 and Nominal and obs based attitude.... 1535899 on one drive
d_p=(y_sc/2)+(l_p/2);

rp1=[0;d_p;0];
rp2=[0;-d_p;0];
Ap1=l_p*b_p;
Ap2=Ap1;

n_p=[0;0;1]; % unit vector that represents the normal to the panel

% Inertia matrix after the panels are deployed

m_p1=3.8*Ap1; % mass of the panels computed as 3.8 kg/m^2 [control it]
m_p2=3.8*Ap2;

Ip1_x=(1/12)*m_p1*(l_p^2) + m_p1*d_p^2;
Ip1_y=(1/12)*m_p1*(b_p^2) + m_p1*d_p^2;
Ip1_z=(1/12)*m_p1*(l_p^2 + b_p^2) + m_p1*d_p^2;

Ip2_x=(1/12)*m_p2*(l_p^2) + m_p2*d_p^2;
Ip2_y=(1/12)*m_p2*(b_p^2) + m_p2*d_p^2;
Ip2_z=(1/12)*m_p2*(l_p^2 + b_p^2) + m_p2*d_p^2;

ip_1=[Ip1_x;Ip1_y;Ip1_z];
ip_2=[Ip2_x;Ip2_y;Ip2_z];

Ip_1_transp=diag(ip_1);
Ip_2_transp=diag(ip_2);

I_tot=I+Ip_1_transp+Ip_2_transp;

I_tot_inv=I_tot\eye(3,3);

I_max = max(diag(I_tot));
I_min = min(diag(I_tot));

%% Orbit Data:

Re = astroConstants(23); %km
mu = astroConstants(13); % [km^3 / s^2]
a = 7714.4278; % [km] AVISO
e = 0.000095; % [-] AVISO
R = a; % [km] Orbit Radius equal to a-> circular approx
h = R-Re; % [km] altitude of orbit
n = sqrt(mu/ a^3); % mean angular velocity (rad/s)
T = (2*pi)/n; % Period (s)
we=15.04*(2*pi/360); %[rad/h] ang vel of Earth
we=we/3600; % rad/s ang vel of Earth
v_circ = sqrt(mu/R); % [km/s] speed of orbit, circular


F_s = 1358 + 500 + 117; %[W/m^2] Direct solar radiation + Radiation refl from E + Albedo at 1000 km slides di Bernelli

Cg_mis = 0.03; %[m] Slides lavagna, misalignment of centre of mass [worst case]

q_p = 0.6; % reflectivity coeff [-] assumed from slides, worst case 0.5-0.6

rho_airh = 3.019e-15; % [kg/m^3] density of air at our altitude of 1336 km (model stops at 1000, worst case assumed at 1000) [check]

Cd = 2.5; % [-] worst case 2-2.5 [check altitude]

v_r = 1000*(v_circ - we*R); % [m/s] relative velocity

D = (1.4*10^-3)*Mass_sc; %[A*m^2] Assume Class II and spinning, magnetic torques comperable to other torques [Bernelli]
% D is 1.4 for spinning, 3.5 for not spinning. Lavagna said D between 1 - 20
M = 7.96e15; %[Tm^3] Assumed from slides

B = 2*M/((R*1000)^3); % [T] 

%% Disturbances:

T_gg = ((3*astroConstants(13))/(2*R^3))*(I_max - I_min); % [Nm] [assumed theta 45 degrees, so worst case]

T_SRP = (F_s/(astroConstants(5)*1000))*(Ap1+Ap2+A3)*(1+q_p)*Cg_mis; % [Nm] [hypotezied that you have the area of the 2 panels and of the head]

T_air_drag = 0.5*rho_airh*Cd*A1*(v_r^2)*Cg_mis; % [Nm]

T_mag = D*B; % [Nm]

T_tot = T_gg + T_SRP + T_air_drag + T_mag; %[Nm]

T_tot_real =2*T_tot; % [Nm] 100% margin is applied

%% RW:

% Data:

H_RW = 12; % [Nms] Taken from Reaction Wheel survey, less advenced version of Jason-2 RW
T_RW = 0.075; % [Nm] Taken from Reaction Wheel survey, less advenced version of Jason-2 RW

H_max = T_tot_real*T; %[Nms] Maximum angular momentum stored for 1 orbit period -> max ang. mom. stored for RW

n_orbit_sat = H_RW/H_max;  % [-] number of orbits before the wheels are saturated

theta_slew_max = pi; %[rad] Assumed worst case for a slew

slew_rate = deg2rad(0.5); %[rad/s] usual value, if nothing else is known

t_slew_man = theta_slew_max/slew_rate; % [s] Assumed almost immediate slew

T_max_slew = 4*theta_slew_max*I_max/(t_slew_man^2); % [Nm] 


% if the T required by the slew is higher then recompute the t and
% theta_dot to have a feasible slew with those RW

if T_max_slew > T_RW

    t_slew_man_real = sqrt(4*theta_slew_max*I_max/T_RW); % [s]

    slew_rate_real = rad2deg(theta_slew_max/ t_slew_man_real); % [deg/s]

    H_slew = T_RW*t_slew_man_real; % [Nms]

else 
   
    t_slew_man_real = t_slew_man; % [s]
    slew_rate_real = rad2deg(theta_slew_max/ t_slew_man_real);  % [deg/s]

    H_slew = T_max_slew*t_slew_man_real; % [Nms]

end

D_desat_slew = 180; %[Am^2]

t_desat_slew = H_slew/(B*D_desat_slew); %[s]

%% Desaturation of Disturbances with Magnetic Torquers:

D_magt = 60; % [A*m^2] 1996 data paper of Jason

n_mag = 3; %[-] number of magneto torquers: 3 available, 1 for conservative

t_desat_magt = H_RW/(D_magt*B); % [s] 

t_desat_perc = t_desat_magt/T; % [-]

% n_s = 1e5; %[-] number of spires, GUESSED VALUE, PLS FIND IT
% 
% S =(pi*0.1^2)/4; % [m^2] area of coil, hyp: 10 cm of diam  GUESSED VALUE, PLS FIND IT
% 
% I_desat = D_magt/(S*n_s); %[A] current

%% Slew with Magneto torquers

T_mag_torq = 20*B; % [Nm]

theta_point = pi; %[rad] 0.1 deg from Pointing budget

theta_point_rate = deg2rad(0.5); %[rad/s] hyp, cannot go too fast

t_point_man = theta_point/theta_point_rate; % [s] Assumed almost immediate slew

T_max_point = 4*theta_point*I_max/(t_point_man^2); % [Nm] 

if T_max_point > T_mag_torq

    t_point_man_real = sqrt(4*theta_point*I_max/T_mag_torq); % [s]

    theta_point_rate_real = rad2deg(theta_point/ t_point_man_real); % [deg/s]

    H_point = T_mag_torq*t_point_man_real; % [Nms]

else 
   
    t_point_man_real = t_point_man; % [s]
    theta_point_rate_real = rad2deg(theta_point/ t_point_man_real);  % [deg/s]
    H_point = T_max_point*t_point_man_real; % [Nms]
end

% D_desat_point = 20; %[Am^2]
% 
% t_desat_point = H_point/(D_desat_point*B); % [s]


% %% Desaturation with Thrusters:
% 
% T_thr = 1; % [N] Thrust provided by 1 thruster
% 
% L = 0.5; % [m] hyp: lever arm of thrusters
% 
% n_thrusters = 3; %[-]
% 
% t_desat_T = H_max/(n_thrusters*L*T_thr); %[s]
% 
% Isp = 220; %[s]
% 
% m_p_desat = (t_desat_T*T_thr)/(9.81*Isp); % [kg]
% 
% m_p_desat_tot = m_p_desat*5*365*24*3600/T; % [kg] total m prop to be used if desaturation with thrusters (nominal mission: 5y)




