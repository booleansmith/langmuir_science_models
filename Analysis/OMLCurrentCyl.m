function [I] = OMLCurrentCyl(x,phi)
% Function NAME:
%
%   OML Current Cylindrical
% 
% DESCRIPTION:
%
%   Produce a current estimation based on an applied potential, phi, and
%   plasma parameters, defined by X. The OML current approximation is based
%   on Swenson's notes, "Lecture 12. OML Collection Currents" 
%   Eqnuations 12.58-12.60. This function is intended to work with
%   LSQ curve fitting to determine plasma paramaters Temperature, and
%   density
%
% INPUTS:
%
%   x   - (double array) a vector of plasma paramters
%       x(1): n      - number density (m^-3)
%       x(2): T      - temperature (K)
%       x(3): phi_p - plasma potential (not used)
%       x(4): Ap    - Surface Area of the probe (m^2)
%       x(5): Rp    - Radius of the probe (m)
%
%   phi - (double array) a vector of applied bias potentials. This should
%       represent a best guess for the potential across the sheath
%
%
% OUTPUTS:
%
%   I - (double array) the current collected for each applied bias
%       potential and plasma parameters
%   n - (double) the density estimate produced by the function
% 
% ASSUMPTIONS AND LIMITATIONS:
%
%   none
%
% Jason Powell
% Feb 01, 2023

%% Constants
k_b = 1.38e-23;  % Boltsmann constant m^2*kg*s^-2*K^-1
epsilon_0 = 8.8541872e-12; %Permittivity of free space m^-3kg^-1s^4A^2
e   = 1.602e-19; % elementary charge
beta = 0.5;

%% Calculate the Debye Length and probe ratio, epsilon

lambda_debye = sqrt(epsilon_0*k_b*x(2)/(x(1)*e^2));
epsilon = (lambda_debye + x(5))/x(5);

%% Oxygen Ion current collection

m_i = 2.65e-26; % O+ ion mass, kg

Jsat_i = e*x(1)*sqrt(k_b*x(2)/(2*pi*m_i));
Phi_si = -e*phi/(k_b*x(2));
Phi_si_tilde = sqrt(Phi_si/(epsilon^2-1));

F_i = zeros(length(phi),1);
F_i(Phi_si < 0) = exp(Phi_si(Phi_si < 0));
F_i(Phi_si >= 0) = epsilon*erf(Phi_si_tilde(Phi_si >= 0)) ... 
                   +exp(Phi_si(Phi_si >= 0)).*erfc(epsilon*Phi_si_tilde(Phi_si >= 0));

I_i = Jsat_i*x(4)*F_i;

%% Electron Current

m_e = 9.109e-31; % electron mass, kg

Jsat_e = -e*x(1)*sqrt(k_b*x(2)/(2*pi*m_e));
Phi_se = e*phi/(k_b*x(2));
Phi_se_tilde = sqrt(Phi_se/(epsilon^2-1));

F_e = zeros(length(phi),1);
F_e(Phi_se < 0) = exp(Phi_se(Phi_se < 0));
F_e(Phi_se >= 0) = epsilon*erf(Phi_se_tilde(Phi_se >= 0)) ... 
                   +exp(Phi_se(Phi_se >= 0)).*erfc(epsilon*Phi_se_tilde(Phi_se >= 0));

I_e = Jsat_e*x(4)*F_e;
    
%% Sum

I = -I_e - I_i;





end