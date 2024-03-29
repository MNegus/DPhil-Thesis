function SubstrateFunctions = platesubstratefunctions(ts, ws, w_ts, w_tts, epsilon)
%PLATESUBSTRATEFUNCTIONS Defines substrate functions given discrete
%solution 
%   Given a discrete solution in terms of the time (ts) and the arrays for
%   the plate solution (ws, w_ts, w_tts), we define a SubstrateFunctions
%   struct that has interpolated functions for w and the substrate
%   dependent quantities such as the turnover point. 

    %% Set parameters
    SubstrateFunctions.dimension = "axi";
    SubstrateFunctions.epsilon = epsilon;
    
    %% Define interpolated functions for rigid displacement term
    SubstrateFunctions.a = @(t) interp1(ts, ws, t, 'linear', 'extrap');
    SubstrateFunctions.a_t = @(t) interp1(ts, w_ts, t, 'linear', 'extrap');
    SubstrateFunctions.a_tt = @(t) interp1(ts, w_tts, t, 'linear', 'extrap');
    
    % Curved displacement terms set to 0
    SubstrateFunctions.b = @(t) zeros(size(t));
    SubstrateFunctions.b_t = @(t) zeros(size(t));
    SubstrateFunctions.b_tt = @(t) zeros(size(t));
    
    % Set w functions (identical to a, fix this later)
    SubstrateFunctions.w = SubstrateFunctions.a;
    SubstrateFunctions.w_t = SubstrateFunctions.a_t;
    SubstrateFunctions.w_tt = SubstrateFunctions.a_tt;
    
    %% Load substrate dependents
    SubstrateFunctions = substratedependents(SubstrateFunctions);

end