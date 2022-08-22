function Fs = innerforce(ts, SubstrateFunctions)
%%innerforce
%   Returns the inner solution for the force on the substrate in the
%   two-dimensional, plate impact case.

    % Load in functions
    epsilon = SubstrateFunctions.epsilon;
    ds = SubstrateFunctions.d(ts);
    d_ts = SubstrateFunctions.d_t(ts);
    Js = SubstrateFunctions.J(ts);

    % Initialise Fs
    Fs = zeros(size(ts));
    
    % Extract quantities from t > 0
    ds = ds(ts > 0);
    d_ts = d_ts(ts > 0);
    Js = Js(ts > 0);

    % Solve for c(t)
    options = optimoptions('fsolve', 'OptimalityTolerance', 1e-12, ...
        'StepTolerance', 1e-12, 'Display', 'off');
    cGuess = 2 * sqrt(pi * ds ./ (4 * epsilon^2 * Js));
    zero_fun = @(c) c.^2 + 4 * c + 2 * log(c) + 1 ...
        - pi * ds ./ (epsilon^2 * Js);
    cs = fsolve(zero_fun, cGuess, options);
    
    %% Return forces
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional force
        
        % At t == 0, we have that the inner force is equal to the overlap
        zeroIdx = find(ts == 0);
        if ~isempty(zeroIdx)
            Fs(zeroIdx) = 2 * sqrt(2) * SubstrateFunctions.C(0);
        end
        
        Fs(ts > 0) = (8 * epsilon * d_ts.^2 .* Js .* cs) / pi;
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric force
        Fs(ts > 0) = (8 * epsilon^4 * d_ts.^2 .* Js.^2 .* cs / pi) ...
            .* (pi * ds ./ (epsilon^2 * Js) - cs.^2 / 3 - 2 * cs - 2 * log(cs) +1);
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end

end