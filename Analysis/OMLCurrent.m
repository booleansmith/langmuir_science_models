function [I] = OMLCurrent(x,phi)
% Function NAME:
%
%   OML Current
% 
% DESCRIPTION:
%
%   Produce a current estimation based on an applied potential, phi, and
%   plasma parameters, defined by X. This function is intended to work with
%   LSQ curve fitting to determine plasma paramaters Temperature, and
%   density
%
% INPUTS:
%
%   x   - (double array) a vector of plasma paramters
%       x(1): n - number density
%       x(2): T - temperature
%       x(3): phi_p - plasma potential (not used)
%       x(4): beta - OML parameter
%       x(5): Ap- Area of the probe
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
e   = 1.602e-19; % elementary charge


%% Oxygen Ion current collection

m_i = 2.65e-26; % O+ ion mass, kg

Jsat_i = e*x(1)*sqrt(k_b*x(2)/(2*pi*m_i));
Phi_si = -e*phi/(k_b*x(2));

F_i = zeros(length(phi),1);
F_i(Phi_si < 0) = exp(Phi_si(Phi_si < 0));
F_i(Phi_si >= 0) = (1+Phi_si(Phi_si >= 0)).^x(4);

I_i = Jsat_i*x(5)*F_i;

%% Electron Current

m_e = 9.109e-31; % electron mass, kg

Jsat_e = -e*x(1)*sqrt(k_b*x(2)/(2*pi*m_e));
Phi_se = e*phi/(k_b*x(2));

F_e = zeros(length(phi),1);
F_e(Phi_se < 0) = exp(Phi_se(Phi_se < 0));
F_e(Phi_se >= 0) = (1+Phi_se(Phi_se >= 0)).^x(4);

I_e = Jsat_e*x(5)*F_e;
    
%% Sum

I = -I_e - I_i;





end