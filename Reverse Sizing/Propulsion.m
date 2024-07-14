close all; clear all; clc; 
% SSEO PROJECT - SIZING PROPULSION

% REQUIREMENTS

% DV requested
total_DV=240; % 130 m/s with 100% margin (stochastic maneuvers)

List_ISP=[50,220,310]; % [s] list of specific impulse per propellant type 
% [cold gas, monoprop,biprop]
% We don't use solid fuel because we want throttling, and reignition
g0=9.80625; %[m/s^2] gravity acceleration
List_MR=exp(total_DV./(List_ISP*g0));
MR_mono=List_MR(2);
% We do not take cold gas because need of throtting for precise maneuvers
% MR for monopropellant and bipropellant are similar
% but monopropellant solution is much more simple than bi prop
% and thrust needed is not a matter (explain why)

% Computation of masses and propellant properties
m_pl=120; % [kg] mass of Payload (from Handbook)
m_dry=420; % [kg] dry mass at launch for an Earth-orbiting s/c
m_dry_real=m_dry*1.2;% [kg] real dry mass with 20% margin
m_prop=m_dry_real*(MR_mono-1); % [kg] propellant mass
m_prop_real=m_prop*1.055; % [kg] real propellant mass (2% margin for ullage, 3% for residuals, 0.5% for uncertainty)
rho_prop=1.01*1e3; % [kg/m^3] density of hydrazine
V_prop=m_prop_real/rho_prop; %[m^3] volume of propellant needed
V_prop_real=V_prop*1.1; %[m^3] volume of propellant with 10% margin for unusable volume;

% Define B the blowndown ratio as the unknown
nbr_points=100;

% Creating Lists
List_B=linspace(1,6,nbr_points); % [-] list of blowdown number contained between 4-6 
List_V_gas_i=zeros(nbr_points,1);
List_V_gas_f=zeros(nbr_points,1);
List_V_tank=zeros(nbr_points,1);
List_P_gas_i=zeros(nbr_points,2);
List_P_gas_f=zeros(nbr_points,2);
List_m_press_He=zeros(nbr_points,2);
List_m_press_N2=zeros(nbr_points,2);
List_r_tank=zeros(nbr_points,1);
List_t_tank=zeros(nbr_points,1);
List_m_tank=zeros(nbr_points,1);
List_dP_feed=zeros(nbr_points,2);

% Defining the pressure system characteristics
dP_feed=50*1e3*[1,1]; % [Pa] 
P_inlet=[5.5,22]*1e5; %[Pa] feed pressure of the thruster
%P_inj=0.3*P_chamber; %[Pa] injector delta pressure
R_spec_He=2077.3; %[J/kg/K] specific perfect gas constant for He
R_spec_N2=296.8; %[J/kg/K] specific perfect gas constant for N2
T_tank=293; %[K] Temperature inside the tank
P_tank_i=22*1e5*[1,1];% [Pa] initial pressure of propellant
Delta_pressure=dP_feed+P_inlet; %[Pa] minimum pressure in tank  


% Properties of the tank
sigma=880*1e6; % strength Ti6A14V
rho_tank=4420; % density Ti6A14V


n_case=1;% 1 is best case, 2 worst case
% Start a counter
n=0;
for i = 1:nbr_points
    B=List_B(i);
    % By considering isothermal conditions
    V_gas_i= V_prop_real/(B-1);
    V_gas_f= V_gas_i*B;
    P_tank_f=[0,0];
    if P_tank_i(n_case)/B<Delta_pressure(n_case)
        return
    end
    P_tank_f=P_tank_i/B;
    n=n+1;
    m_press_He=P_tank_i*V_gas_i/(R_spec_He*T_tank);
    m_press_N2=P_tank_i*V_gas_i/(R_spec_N2*T_tank);
    V_tank=V_gas_i+V_prop_real;
    r_tank=(V_tank/(4*pi/3))^(1/3);
    t_tank=58.2*1e5*1.5*r_tank/(2*sigma);% PROBLEM WITH THE thickness
    m_tank=rho_tank*4/3*pi*((r_tank+t_tank)^3-r_tank^3); % PROBLEM WITH THE MASS
    dP_feed=P_tank_f-P_inlet;
    List_V_gas_i(i)=V_gas_i;
    List_V_gas_f(i)=V_gas_f;
    List_P_gas_f(i,:)=P_tank_f;
    List_P_gas_i(i,:)=P_tank_i;
    List_m_press_He(i,:)=m_press_He*1.2; % 20% margin for pressure-fed systems
    List_m_press_N2(i,:)=m_press_N2*1.2; % 20% margin for pressure-fed systems
    List_V_tank(i)=V_tank*1.01; % 1% margin to consider bladder volume 
    List_r_tank(i)=r_tank;
    List_t_tank(i)=t_tank;
    List_m_tank(i)=m_tank;
    List_dP_feed(i,:)=dP_feed;
end

%% Study for B Closest to 4 so that final pressure in tank is different to 0
clc;
fprintf("For B = " + string(List_B(n)) + "\n")
fprintf('The volume of the tank [L]  is: \n')
List_V_tank(n)*1e3
fprintf('The initial volume of the pressurant gas [L]  is: \n')
List_V_gas_i(n)*1e3
fprintf('The final volume of the pressurant gas [L]  is: \n')
List_V_gas_f(n)*1e3
fprintf('The initial pressure in the tank [bar]  is: \n')
List_P_gas_i(n,n_case)*1e-5
fprintf('The final pressure in the tank [bar]  is: \n')
List_P_gas_f(n,n_case)*1e-5
fprintf('The mass of the pressurant [kg] (for He) \n')
List_m_press_He(n,n_case)
fprintf('The mass of the pressurant [kg] (for N2) \n')
List_m_press_N2(n,n_case)
fprintf('The mass of propellant [kg]: \n')
m_prop_real

% By considering a SPHERICAL tank with Ti6A14V
fprintf('The diameter of the tank [mm]: \n')
List_r_tank(n)*2*1e3
fprintf('The thickness of the tank [mm]: \n')
List_t_tank(n)*1e3
fprintf('The mass of the tank [kg]: \n')
List_m_tank(n)

m_thruster=0.29; %[kg]mass of thruster with valve
m_PS=List_m_tank(n)+List_m_press_N2(n)+4*m_thruster;
fprintf('The mass of the propulsion system [kg]: \n')
m_PS*1.1 % with 10% margins to account for cables etc
fprintf('The power required of the propulsion system [W]: \n')
power_PS=4*(6.4+6.5+2*9.5)

% Feeding line 
mu=0.876*1e-3;%[PA.s] dynamic viscosity of hydrazine https://apps.dtic.mil/sti/tr/pdf/AD0042818.pdf
d=6.25*1e-3; % [m] it is a 1/4"=6.25 mm 
A=pi/4*(d)^2; % [m^2] cross-section area of the pipe line
v_pipe=sqrt(List_dP_feed(n,n_case)*2/rho_prop);% v_pipe=m_pipe/(A*rho_prop);
fprintf('The flow velocity in the pipe [m/s]: \n')
v_pipe
fprintf('The Reynolds number of the  flow [-]: \n')
Re=rho_prop*v_pipe*d/mu
fprintf('The friction factor [-]: \n')
f=colebrook(Re) % friction factor for laminar flow
fprintf('Length of the pipe [m]: \n')
L=List_dP_feed(n,n_case)*2*d/(f*rho_prop*v_pipe^2)

%% Plots 
close all; clc; 
figure;
%plot(List_B,List_m_press_N2(:,2))

%Initial pressure in tank
subplot(2,3,1);
plot(List_B,List_P_gas_i(:,n_case),'DisplayName',"Initial pressure in tank")
set(gca,'fontsize',12)
grid on
legend show

%Final pressure in tank
subplot(2,3,2);
plot(List_B,List_P_gas_f(:,n_case),'DisplayName',"Final pressure in tank")
set(gca,'fontsize',12)
grid on
legend show
%
% plot(List_B,List_m_press(:,1))
% plot(List_B,List_m_press(:,2))

%Volume of tank
subplot(2,3,3);
plot(List_B,List_V_tank*1e3,'DisplayName',"Volume of tank")
set(gca,'fontsize',12)
grid on
legend show

%  Diameter of tank
subplot(2,3,4);
plot(List_B,List_r_tank*2*1e3,'DisplayName',"Diameter of tank")
set(gca,'fontsize',12)
grid on
legend show

%  FLow velocity
subplot(2,3,5);
plot(List_B,sqrt(List_dP_feed(:,n_case)*2/rho_prop),'DisplayName',"Flow velocity in pipe")
set(gca,'fontsize',12)
grid on
legend show