
addpath(fullfile('..','utilities'))
boxpath = getboxpath;
load(fullfile(boxpath,"LTspice","langmuir","SPORTmodelV5.1.mat"),"data");
load(fullfile(boxpath,"Nascap","Data","nascapdata.mat"),"nascapTable");


%% Constants

q_e = -1.602e-19;
m_e = 9.109e-31;
k_b = 1.38e-23;

%% Plot Matlab Data

areas = unique(data.Curarea);
densities = unique(data.Curdensity);
temps = unique(data.Curtemp);

for i = 1:height(areas)
    T0 = data(data.Curarea == areas(i),:);
    figure
    tiledlayout('flow')
    count = 0;
    for j = 1:height(temps)
        
        T1 = T0(T0.Curtemp == temps(j),:);
        for k = 1:height(densities)
            T2 = T1(T1.Curdensity == densities(k),:);
            count = count + 1;
            ax(count) = nexttile;
            plot(T2.v1,T2.("Ix(LP:Surface)"))
            yyaxis right
            plot(T2.v1,T2.("V(lp)"))
            title(sprintf('CurTemp = %0.2f, CurDensity= %0.2e',temps(j),densities(k)))
        end
    end
    sgtitle(sprintf('CurArea = %0.2f', areas(i)))
%     linkaxes(ax,'xy')
end

%% Plot Nascap-2k Data
load(fullfile(boxpath,"Nascap","Data","nascapdata.mat"),"nascapTable");
nascapTable.Conductor_7_Current = nascapTable.Conductor_7_Current*352e-6;

densities = unique(nascapTable.density);
temps = unique(nascapTable.temperature);


figure;
tiledlayout('flow')

count = 0;
for j = 1:height(temps)
    
    T1 = nascapTable(nascapTable.temperature == temps(j),:);
    for k = 1:height(densities)
        T2 = T1(T1.density == densities(k),:);
        T2 = T2(3:end,:);
        count = count + 1;
        % Plot the Current
        ax1(count) = nexttile;
        plot(T2.Conductor_7_Potential - T2.Conductor_1_Potential,T2.Conductor_7_Potential)
        
        title(sprintf('CurTemp = %0.2f, CurDensity= %0.2e',temps(j),densities(k)))
        ylabel('Sheath Potential')

        % Plot the Potentials
        yyaxis right
        plot(T2.Conductor_7_Potential - T2.Conductor_1_Potential,-T2.Conductor_7_Current)

        legend('Potential','Current')
        xlabel('Bias Potential')
        ylabel('Current Into Probe')
        
    end
end
sgtitle('Nascap-2k Currents')
% linkaxes(ax1)


%% Plot together

load(fullfile(boxpath,"LTspice","langmuir","SPORTmodelV5.1.mat"),"data");
load(fullfile(boxpath,"Nascap","Data","nascapdata.mat"),"nascapTable");
nascapTable.Conductor_7_Current = nascapTable.Conductor_7_Current*342e-6;

areas = unique(data.Curarea);
densities = unique(data.Curdensity);
temps = unique(data.Curtemp);


for i = 1:height(areas)
    T0 = data(data.Curarea == areas(i),:);
    figure
    tiledlayout('flow')
    count = 0;
    for j = 1:height(temps)
        NT1 = nascapTable(nascapTable.temperature == temps(j),:);
        T1 = T0(T0.Curtemp == temps(j),:);
        for k = 1:height(densities)
            jsat = -q_e*densities(k)*sqrt(k_b*temps(j)/(2*pi*m_e))*171e-6;
            NT2 = NT1(NT1.density == densities(k),:);
            T2 = T1(T1.Curdensity == densities(k),:);
            NT2 = NT2(3:end,:);
            count = count + 1;

            ax(count) = nexttile;
            semilogy(T2.v1,abs(T2.("Ix(LP:Surface)")))
            hold on
            semilogy(NT2.Conductor_7_Potential -NT2.Conductor_1_Potential,abs(-NT2.Conductor_7_Current))
            yline(jsat)
            title(sprintf('CurTemp = %0.2f, CurDensity= %0.2e',temps(j),densities(k)))
            legend('LTSpice','Nascap-2k','Location','southeast')
        end
    end
    sgtitle(sprintf('CurArea = %0.2f', areas(i)))
%     linkaxes(ax,'xy')
end








