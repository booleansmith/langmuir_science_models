% PROCESS NAME:
%
%   Langmuir Probe Data
% 
% DESCRIPTION:
%
%   Produce syntetic Langmuir probe data to analyse
%
% INPUTS:
%
%   none
%
% OUTPUTS:
%
%   data - (table) contains entries for LP current measurements and 
%          respective voltage potential levels
% 
% ASSUMPTIONS AND LIMITATIONS:
%
%   none
%


% Jason Powell 09 September 2022

% The first model we consider is based on Debchoudhury's method for looking
% at FPMU data on the ISS: "Observations and Validation of Plasma Density,
% Temperature, and O+ Abundance From a Langmuir Probe Onboard the ISS"

% I = -I_ram_i - I_OML_i + I_e
%   I_ram_i = n_j*A_proj*e*V_sc ; The ions that run into collector body
%   I_OML_i = I_th_j*({1-e(phi-phi_p)}/k_b*T_i)^beta ; the OML ion current
%   I_e     = I_th_e*exp({e(phi-phi_p)}/k_b*T_e) ; the electron current 
%   I_th_j  = n_j*q_j*A*sqrt(k_b*T_j/(2*pi*m_j)) ; the thermal current for ions/electrons

% KNOWN CONSTANTS
k_b = 1.38e-23;  % Boltsmann constant m^2*kg*s^-2*K^-1
e   = -1.602e-19; % electron charge coulombs
m_e = 5.486e-4; % electron mass, kg

% Environment Parameters
n = 1e12; % density m^-3
T = 1160; % Temperature, K
m_i = 16; % O+ ion mass

phi = -2:0.01:3; % Bias potential, Volts
phi_p = -1; % plasma potential, Volts
V_sc  = -7.4e3; % SC velocity, m/s


% Probe parameters
A_proj = 1e-3*85.2e-3; % projected area, m^2
A      = pi*0.5e-3*85.2e-3; % OML collection area, m^2
beta = 1;

I_ram = n*A_proj*e*V_sc;
I_th_e = n*e*A*sqrt((k_b*T)/(2*pi*m_e));
I_th_i = -n*e*A*sqrt((k_b*T)/(2*pi*m_i));
I_OML_i = I_th_i*(1-((e*(phi-phi_p))/(k_b*T))).^beta;
I_e     = I_th_e*exp((e*(phi-phi_p))/(k_b*T));

I = -I_ram -I_OML_i + I_e;

figure
plot(phi,I)
yyaxis right
plot(phi(2:end),diff(I))
























