


Phi_s = -3:0.001:5; % Normalized sheath potential (unitless)
epsilon = 1; % collection factor (sheath distance/ radius of the probe)




f = figure;
% tiledlayout(1,2)
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
% nexttile
h1 = stepplot(Phi_s,cylindrical(Phi_s,epsilon));

% nexttile
% h2 = stepplot(Phi_s,spherical(Phi_s,epsilon));




% UI Control
b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',epsilon, 'min',0, 'max',10000);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','10000','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','Epsilon','BackgroundColor',bgcolor);

b.Callback = @(es,ed) updateSystem(h1,spherical(Phi,es.value)); 

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