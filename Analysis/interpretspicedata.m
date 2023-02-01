



addpath(fullfile('..','utilities'))
boxpath = getboxpath;

%% Langmuir Temp

%% Load temperatureData from .txt file, 
% This only needs to be done once, but is nice to see what was done to get the spice data as a table

temperatureData = fullfile(boxpath,"LTspice","langmuir","LangmuirTemp.txt");
data = loadspicefile(temperatureData);
data.appliedpotential = data.v1 + data.("V(scground)");
save(fullfile(boxpath,"LTspice","langmuir","LangmuirTemp.m"),"data")

%% Plot the data

temps = unique(data.steps);

for i = 1:height(temps)
    d2 = data(data.steps == temps(i),:);
    figure
    tiledlayout(2,1)
    % current vs bias voltage
    nexttile
    plot(d2.v1,d2.("Ix(LP:Plasma)"))
    title('Current Voltage Curve')
    xlabel('Bias Voltage [V]')
    ylabel('Current [A]')

    % Applied potential vs bias potential
    nexttile
    plot(d2.v1,d2.v1 + d2.appliedpotential)
    title('Bias Potential vs Applied Potential')
    xlabel('Bias Potential [V]')
    ylabel('Potential Applied to Plasma [V]')

    sgtitle(sprintf('Plasma Temperature: %2.1fK',d2.steps(1)))

end







%% SPORT Model








