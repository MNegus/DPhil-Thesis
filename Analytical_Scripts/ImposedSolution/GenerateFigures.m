clear;
close all;

%% Create figure directories
for dirName = ["Two-dimensional_Figures", "Axisymmetric_Figures"]
    mkdir(dirName);
    mkdir(dirName + "/png");
    mkdir(dirName + "/eps");
    mkdir(dirName + "/fig");
end

%% Loop over dimensions
% for dimension = ["2D", "axi"]
for dimension = "axi"

    %% Time dependents
    TimeDependentsPlot(dimension);
    
    %% Jet evolution
%     JetEvolutionPlot(dimension);
    
    %% Composite pressure 
%     CompositePressurePlot(dimension);
    
    %% Pressure evolution
%     PressureEvolutionPlot(dimension);

    %% Outer free surface plot
%     OuterFreeSurfacePlot(dimension);

    %% Outer outer pressure
%     OuterOuterPressurePlot(dimension);


end

%% Single dimension plots
% OuterStreamlinePressure2D();