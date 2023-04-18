


addpath(fullfile('..','utilities'))
boxpath = getboxpath;
fpath = fullfile(boxpath,"Nascap","Data");

%% 

runs = struct2table(dir(fullfile(fpath,'*.txt')));
runs = runs(contains(runs.name,'Potentials'),:);
run_names = string(runs.name);


for i = 1:height(run_names)
    potentialName = runs.name{i};
    currentName = strrep(runs.name{i},'Potentials','Currents');
    PTable = readtable(fullfile(runs.folder{i},potentialName));
    PTable.Properties.VariableNames(2:end) = strcat(PTable.Properties.VariableNames(2:end),"_Potential");
    CTable = readtable(fullfile(runs.folder{i},currentName));
    CTable.Properties.VariableNames(2:end) = strcat(CTable.Properties.VariableNames(2:end),"_Current");

    T = join(PTable,CTable);


    nameparts = strsplit(run_names(i,1),{'_','.'});
    curN = str2double(nameparts(2));
    curT = str2double(nameparts(3));


    T.density = ones(height(T),1)*curN;
    T.temperature = ones(height(T),1)*curT;

    if i == 1
        nascapTable = T;
    else
        nascapTable = [nascapTable;T];
    end
  
end

save(fullfile(fpath,'nascapdata.mat'),"nascapTable");

