


Phi_s = -3:0.001:5; % Normalized sheath potential (unitless)
epsilon = 1; % collection factor (sheath distance/ radius of the probe)




f = figure;
% tiledlayout(1,2)
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
sys = cylindrical(Phi_s,epsilon);
% nexttile
h1 = stepplot(ax,sys);

% nexttile
% h2 = stepplot(Phi_s,spherical(Phi_s,epsilon));


 

% Equation for spherical probes
function F_sp = spherical(Phi,epsilon) 
    Phi_t = sqrt(Phi./(epsilon.^2-1)); % phi tilde
    F_sp(Phi <= 0) = exp(Phi(Phi <= 0));
    F_sp(Phi > 0) = epsilon.^2*(1-exp(-1*Phi_t(Phi > 0))) + ...
        exp(-1*Phi_t(Phi > 0));
end


% Equation for cylindrical probes
function F_c = cylindrical(Phi,epsilon) 
    Phi_t = sqrt(Phi./(epsilon.^2-1)); % phi tilde
    F_c(Phi <= 0) = exp(Phi(Phi <= 0));
    F_c(Phi > 0) = epsilon*erf(Phi_t(Phi > 0)) + ...
        exp(Phi(Phi > 0)).*erfc(epsilon.*Phi_t(Phi > 0));
end