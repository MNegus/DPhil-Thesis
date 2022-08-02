epsilon = 1;
ALPHA = 2;
BETA = 0;
GAMMA = 500;
tMax = 0.625;
forceType = "composite";

%% Solve ODE
[ts, ws, w_ts, w_tts] = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, forceType);

%% Plot solution
plot(ts, ws);

%% TODO
% * Make interp functions for w etc and put them in a SubstrateFunctions
% struct
% * Use the struct to call substratedependents to get functions for d etc.
% * Call the substrate force functions to get the outer and composite
% pressure 