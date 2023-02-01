
%% Constants Used 

constant = constants();
kb = constant.kb; % kevin boltzmann constant, J/k
q_e = -1*constant.e; % electron charge, coulombs
q_i = constant.e; % ion charge, coulombs
m_oi = constant.moi; % oxygen ion mass, kg
m_e = constant.me; % electron mass, kg

%%

t = 1500; % Temperature, degrees Kelvin
n = 1e12; % Density, m^-3 
phi_p = -2:0.1:3; % Applied probe potential, Volts

%% Ion Currents

% Saturation Current

Isat = q_i*n*sqrt(kb*T/(2*pi*mi));
% Ram Current








