function [I] = OMLCurrentCylMoving(x,phi)
% Function NAME:
%
%   OML Current Cylindrical Moving
% 
% DESCRIPTION:
%
%   Produce a current estimation based on an applied potential, phi, and
%   plasma parameters, defined by X. The OML current approximation is based
%   on Swenson's notes, "Lecture 12. OML Collection Currents" 
%   Eqnuations 12.62-12.65. This function is intended to work with
%   LSQ curve fitting to determine plasma paramaters Temperature, and
%   density
%
% INPUTS:
%
%   x   - (double array) a vector of plasma paramters
%       x(1): n      - number density (m^-3)
%       x(2): T      - temperature (K)
%       x(3): Ap    - Surface Area of the probe (m^2)
%       x(4): Rp    - Radius of the probe (m)
%       x(5): w     - speed of the probe relative to the plasma (m/s)
%       x(6): theta - the angle between the probe axis and the velocity
%                     vector (radians)
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
M_i = sqrt(m_i*x(5)^2*sin(x(6))^2/(2*k_b*x(2)));

F_i = zeros(length(phi),1);
% Define mask functions for when different collection factors should be used
% m1 = Phi_si <= 0 & M_i < 0.1;
% m2 = Phi_si > 0 & M_i < 4;
m1 = Phi_si <= 0;
m2 = Phi_si > 0;
F_i(m1) = exp(Phi_si(m1)).* ... 
                   (1 + ...
                   (0.5-Phi_si(m1))*M_i^2 + ...
                   (0.25*Phi_si(m1).^2 + 0.25*Phi_si(m1))*M_i^4 ...
                   );
F_i(m2) = (2/sqrt(pi))*sqrt(Phi_si(m2)) + exp(Phi_si(m2)).*erfc(sqrt(Phi_si(m2))) ...
                  + (M_i^2/2)*((2/sqrt(pi))*sqrt(Phi_si(m2))+(1-2*Phi_si(m2)).*exp(Phi_si(m2)).*erfc(sqrt(Phi_si(m2)))) ... 
                  - (M_i^4/8)*((2/sqrt(pi))*sqrt(Phi_si(m2))+(0.5-2*Phi_si(m2)-2*Phi_si(m2).^2).*exp(Phi_si(m2)).*erfc(sqrt(Phi_si(m2))));

I_i = Jsat_i*x(3)*F_i;

%% Electron Current

m_e = 9.109e-31; % electron mass, kg

Jsat_e = -e*x(1)*sqrt(k_b*x(2)/(2*pi*m_e));
Phi_se = e*phi/(k_b*x(2));
M_e = sqrt(m_e*x(5)^2*sin(x(6))^2/(2*k_b*x(2)));

F_e = zeros(length(phi),1);
% m1 = Phi_se <= 0 & M_e < 0.1;
% m2 = Phi_se > 0 & M_e < 4;
m1 = Phi_se <= 0;
m2 = Phi_se > 0;
F_e(m1) = exp(Phi_se(m1)).* ... 
                   (1 + ...
                   (0.5-Phi_se(m1))*M_e^2 + ...
                   (0.25*Phi_se(m1).^2 + 0.25*Phi_se(m1))*M_e^4 ...
                   );
F_e(m2) = (2/sqrt(pi))*sqrt(Phi_se(m2)) + exp(Phi_se(m2)).*erfc(sqrt(Phi_se(m2))) ...
                  + (M_e^2/2)*((2/sqrt(pi))*sqrt(Phi_se(m2))+(1-2*Phi_se(m2)).*exp(Phi_se(m2)).*erfc(sqrt(Phi_se(m2)))) ... 
                  - (M_e^4/8)*((2/sqrt(pi))*sqrt(Phi_se(m2))+(0.5-2*Phi_se(m2)-2*Phi_se(m2).^2).*exp(Phi_se(m2)).*erfc(sqrt(Phi_se(m2))));

I_e = Jsat_e*x(3)*F_e;
    
%% Sum

I = I_e + I_i;





end