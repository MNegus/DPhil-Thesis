function [ts, ws, w_ts, w_tts] = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, forceType)
%PLATESOLUTION Finds solution for plate displacement
%   Finds solution for the plate displacement, w(t), given by
%   (ALPHA / epsilon^2) * w_tt(t) + BETA * w_t(t) + epsilon^2 GAMMA * w(t)
%   = F(t),
%   where F(t) is the force (either given by the outer or composite
%   solution). 

    %% Add relevant paths
    addpath("../");
    addpath("../Forces");
    
    %% Computational parameters
    t0 = 0; % Initial time (instead of initialising at t = 0)
    RelTol = 1e-4; % Relative tolerance of ode15i
    AbsTol = 1e-5; % Absolute tolerance of ode15i
    MaxStep = 1e-3; % Max timestep for ode15i
    tSpan = [t0 tMax]; % Time span

    
    %% ODE residual function
    function res = FullPlateResFun(t, y, yp, ALPHA, BETA, GAMMA, epsilon, forceType)
    %FULLPLATERESFUN Residual function used for ode15i
    %   Residual function for solving the plate ODE, where w -> y and w_t
    %   -> yp.
    %   TODO: Make a separate function to do the loading below to make this
    %   all a bit neater. 
        % Save substrate functions struct
        SubstrateFunctions.dimension = "axi";
        SubstrateFunctions.epsilon = epsilon;
        SubstrateFunctions.a = @(t) y(1);
        SubstrateFunctions.a_t = @(t) yp(1);
        SubstrateFunctions.a_tt = @(t) yp(2);
        SubstrateFunctions.b = @(t) 0;
        SubstrateFunctions.b_t = @(t) 0;
        SubstrateFunctions.b_tt = @(t) 0;
        
        % Load substrate dependents
        SubstrateFunctions = substratedependents(SubstrateFunctions);
        
        % Find substrate force
        if t == 0
            forceTerm = 0;
            
        else
            if forceType == "outer"
                % Outer force
                forceTerm = outerforce(t, SubstrateFunctions);
            elseif forceType == "composite"
                % Composite force
                [forceTerm, ~, ~] = substrateforce(t, SubstrateFunctions);
            else
                error("Invalid forceType. Needs to be 'outer' or 'composite'.");
            end
            
        end
        % Return residual
        res = [yp(1) - y(2); ...
            ALPHA * yp(2) / epsilon^2 + BETA * yp(1) + epsilon^2 * GAMMA * y(1) ...
            - forceTerm];
        
        
    end

    %% Set initial conditions
    w0 = 0;
    w_t0 = 0;
    w_tt0 = 0;
    y0 = [w0; w_t0];
    yp0 = [w_t0; w_tt0];
    
    %% Solve ODE
    odeOptions = odeset('RelTol', RelTol, 'AbsTol', AbsTol, ...
        'MaxStep', MaxStep, 'Stats', 'on');
    PlateResFun = @(t, y, yp) FullPlateResFun(t, y, yp, ALPHA, BETA, GAMMA, epsilon, forceType);
    [ts, ys] = ode15i(PlateResFun, tSpan, y0, yp0, odeOptions);
    
    %% Extract w and its derivatives
    ws = ys(:, 1);
    w_ts = ys(:, 2);
    w_tts = diff(w_ts) ./ diff(ts);
    
    %% Shortens arrays to make same length
    ts = ts(1 : end - 1);
    ws = ws(1 : end - 1);
    w_ts = w_ts(1 : end - 1);


end