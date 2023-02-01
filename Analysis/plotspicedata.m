function plotspicedata(S)

% Plotting sizes
    lrg = 24;
    med = 18;
    sml = 12;

figure

areas = unique(S.Table.scarea);
myLegend =[];
tiledlayout(3,1)

% Plot the Voltage Sweep of Langmuir Probe

SA_probe = 297e-6;


for i = 1:length(areas)
    T = S.Table(S.Table.scarea == areas(i),:);
    myLegend = [myLegend,sprintf("SC Area = %0.3f",areas(i))];

    Vb = T.("V(lp)") - T.("V(scground)"); % biased voltage
    Vp = T.("V(lp)"); % voltage w reference to plasma
    nexttile(1)
        hold on
        plot(Vb,T.("V(lp)"),'LineWidth',2)
        hold off

    nexttile(2)
        hold on
        plot(Vb,T.("Ix(lp:SURFACE)")*1e6,'LineWidth',2)
        hold off

    nexttile(3)
        if (i == 1)
            semilogy(Vb,abs(T.("Ix(lp:PLASMA)"))*1e6,'LineWidth',2)
        else
            hold on
            semilogy(Vb,abs(T.("Ix(lp:PLASMA)"))*1e6,'LineWidth',2)
            hold off
        end

end

nexttile(1)
    T = S.Table(S.Table.scarea == areas(1),:);
    hold on
    plot(Vb, Vb,'LineWidth',3,'LineStyle','--','Color','r');
    hold off
    ylim([-2.5,3.5])
    xlim([-2,2.99])
%     xlim([T.time(1)*1000+1,T.time(end)*1000+10])
    title('Langmuir Probe Potential','FontSize',med)
    xlabel('Probe to Spacecraft Potential (V)','FontSize',med)
    ylabel('Probe to Plasma (V)','FontSize',med)
    legend([myLegend,'Probe Voltage'],'Location','northwest')

nexttile(2)
    title('Current of Langmuir Probe','FontSize',med)
    xlabel('Probe to Spacecraft Potential (V)','FontSize',med)
    ylabel('Current to Plasma (\muA)','FontSize', med)
    xlim([-2,2.99])
    legend(myLegend,'Location','northwest')

nexttile(3)
    title('Current out of Langmuir Probe (Log|i|)','FontSize',med)
    %xlim([T.time(1)*1000+1,T.time(end)*1000+10])
    %ylim([-15 5])
    xlim([-2,2.99])
    xlabel('Probe to Spacecraft Potential (V)','FontSize',med)
    ylabel('Current to Plasma (\muA)','FontSize', med)
    legend(myLegend,'Location','northwest')
% Compare to NASCAP
% NASCAP FILES ARE LOADED SEPERATLY
% TODO: FIX THIS SO THE SCRIPT RUNS WITHOUT OTHER DATA


nexttile(1)
    T = S.Table(S.Table.scarea == areas(1),:);
    hold on
%     plot(Vb, Vb,'LineWidth',3,'LineStyle','--','Color','r');
        plot(potentials.biasVoltage,potentials.Probe,'LineWidth',2) % FROM NASCAP SIMULATION
        plot(potentials.biasVoltage,potentials.biasVoltage,'LineWidth',2)
    hold off
    ylim([-2.5,3.5])
    xlim([-2,2.99])
%     xlim([T.time(1)*1000+1,T.time(end)*1000+10])
    title('Langmuir Probe Potential','FontSize',med)
    xlabel('Probe to Spacecraft Potential (V)','FontSize',med)
    ylabel('Probe to Plasma (V)','FontSize',med)
    legend([myLegend,'Probe Voltage'],'Location','east')
        legend([myLegend,'Nascap Estimate','Nascap Bias'],'Location','northwest')

nexttile(2)
        hold on
        plot(currents.biasVoltage,-1*currents.Probe*SA_probe*1e6,'LineWidth',2)
        hold off
    title('Current of Langmuir Probe','FontSize',med)
    xlabel('Probe to Spacecraft Potential (V)','FontSize',med)
    ylabel('Current to Plasma (\muA)','FontSize', med)
    xlim([-2,2.99])
    legend(myLegend,'Location','east')
        legend([myLegend,'Nascap Estimate'],'Location','northwest')

nexttile(3)
    title('Current out of Langmuir Probe (Log|i|)','FontSize',med)
        hold on
        semilogy(currents.biasVoltage,abs(-1*currents.Probe)*SA_probe*1e6,'LineWidth',2)
        hold off
    %xlim([T.time(1)*1000+1,T.time(end)*1000+10])
    %ylim([-15 5])
    xlim([-2,2.99])
    xlabel('Probe to Spacecraft Potential (V)','FontSize',med)
    ylabel('Current to Plasma (\muA)','FontSize', med)
    legend(myLegend,'Location','east')
        legend([myLegend,'Nascap Estimate'],'Location','northwest')
end