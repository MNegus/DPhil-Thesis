
%% Create figure directories
for dirName = ["Two-dimensional_Figures", "Axisymmetric_Figures"]
    mkdir(dirName);
    mkdir(dirName + "/png");
    mkdir(dirName + "/eps");
    mkdir(dirName + "/fig");
end

for dimension = ["2D", "axi"]

    %% Time dependents
    TimeDependentsPlot(dimension);
    
    %% Jet evolution
    JetEvolutionPlot(dimension);
end