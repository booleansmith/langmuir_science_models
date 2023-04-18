
addpath(fullfile('..','utilities'))
boxpath = getboxpath;
save(fullfile(boxpath,"LTspice","langmuir","SPORTmodelV5.1.mat"),"data")

%%

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
            title(sprintf('CurTemp = %0.2f, CurDensity= %0.2e',temps(j),densities(k)))
        end
    end
    sgtitle(sprintf('CurArea = %0.2f', areas(i)))
    linkaxes(ax,'xy')
end




