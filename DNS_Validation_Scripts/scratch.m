%% Parameters
% Time parameters used by the simulation
L = 2;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

%% Save imposed solution
epsilon = 1; % Analytical parameter for small time 
a = 0.0025; % Imposed parameter
[ws, w_ts, w_tts] = ImposedPlate(ts_analytical, a);

%% Find pressure
lambda = pi / (2 * L);
p0s = AnalyticalPressure2D(ts_analytical, ws, w_ts, w_tts, lambda, epsilon);

%% Plot pressure
plot(ts_analytical, p0s);