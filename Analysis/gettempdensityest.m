function [t,n] = gettempdensityest(V,I,Tg,densityGuess,dodebug)
% Function NAME:
%
%   Get Temperature and Density Estimation
% 
% DESCRIPTION:
%
%   Takes current voltage current measurements and returns temperature and
%   denisty parameters
%
% INPUTS:
%
%   V  - (double array) a vector of applied potentials
%   I  - (double array) a vector of measured currents
%        corresponding to applied potentials
%   tg - (double) the estimated temperature
%   ng - (double) the estimated density
%
% OPTIONAL INPUTS
%
%   dodebug - (boolean) optionally enable plots for debugging
%
% OUTPUTS:
%
%   t - (double) the temperature estimate produced by function
%   n - (double) the density estimate produced by the function
% 
% ASSUMPTIONS AND LIMITATIONS:
%
%   none
%
% Jason Powell
% January 31, 2023
%%

if ~exist('dodebug','var')
    dodebug = false;
end


% KNOWN CONSTANTS
k_b = 1.38e-23;  % Boltsmann constant m^2*kg*s^-2*K^-1
e   = 1.602e-19; % elementary charge
q_e   = e; % electron charge coulombs 
q_i   = e;  % ion charge coulombs
m_e = 9.109e-31; % electron mass, kg

% Probe parameters
probeRadius = 0.5e-3;    % meters
probeHeight = 85.2e-3; % meters
A_proj = probeRadius*2*probeHeight; % projected area, m^2
A      = pi*probeRadius*probeHeight; % OML collection area, m^2
beta = 0.5;
V_sc  = 7.4e3; % SC velocity, m/s
m_i = 2.65e-26; % O+ ion mass, kg


%% Pre-treat the data a little bit

[I,ia,~] = unique(I,'stable');
V = V(ia);

%% Solve for stage 1 parameters

% floating potential
vf = interp1(I,V,0); % The voltage at the zero current threshold

% plasma potential
Idot = gradient(I,V); % di/dv
phi_p = max(Idot(V > vf & V < vf+1.5));

% ion density
index = V > -2 & V < -1.5;
p = polyfit(V(index),I(index),1);
I_ram_est = -polyval(p,phi_p);
ni_est = I_ram_est/(A_proj*e*V_sc);


if dodebug 
    figure
    plot(V,I)
    hold on
    plot(vf,0,"*",'MarkerSize',5)
    plot(V,polyval(p,V),'--k','LineWidth',1)
    xline(phi_p,'--k')
    plot(phi_p,-I_ram_est,"*",'MarkerSize',5)
    text(phi_p+.05,-I_ram_est-2.5e-7,sprintf('n_i = %0.2e',ni_est))
    ylabel('Return Current (A)')
    yyaxis right
    plot(V,Idot)
    hold on
    plot(phi_p,Idot(Idot == phi_p),"*",'MarkerSize',5)
    xlabel('Bias Potential (V)')
    ylabel('dI/dV')

    title('First Stage parameters from IV curve')
    legend('I','Vf','','','I_{Ram}','dI/dV','\Phi_p')
end


%% Do non-linear curve fitting

% x(1): n - number density
% x(2): T - temperature
% x(3): phi_p - plasma potential
% x(4): beta - OML parameter
fun = @(x,phi) +x(1)*A_proj*e*V_sc ... % ion ram current
               +x(1)*e*A*sqrt((k_b*x(2))/(2*pi*m_i))*(1-((e*phi)/(k_b*x(2)))).^x(4) ... % ion current
               +x(1)*-e*A*sqrt((k_b*x(2))/(2*pi*m_e))*exp((e*phi)/(k_b*x(2))); % electron current
 


% for now we use known values for phi_p and beta
T_guess = 350; % kelvin
phi_guess = phi_p;
% beta_guess = beta;
beta_guess = 1;

% x0 = [ni_est,T_guess,phi_guess,beta_guess]; % initial guesses for lsqcurvefit
% bounds(1,:) = x0 - 0.5*x0;
% bounds(2,:) = x0 + 0.5*x0;
% tolerance = 1e-25;
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','StepTolerance',tolerance,'FunctionTolerance',tolerance);
% x = lsqcurvefit(fun,x0,V,I,bounds(1,:),bounds(2,:),options);
x2 = [1e12,1200,phi_p,0.5,A];

figure
plot(V,I)
hold on
plot(V,OMLCurrent(x2,V))
yyaxis right
plot(V,fun(x2(1:4),V))
legend('Original IV','Swenson Estimate','Debchoudhury estimate')


end

