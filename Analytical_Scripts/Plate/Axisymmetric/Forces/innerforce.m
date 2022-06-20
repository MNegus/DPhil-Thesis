function Fs = innerforce(ts, ds, d_ts, Js, epsilon)
%%innerforce
%   Returns the inner solution for the force on the substrate in the
%   two-dimensional, plate impact case.

    % Initialise Fs
    Fs = zeros(size(ts));
    
    % At t == 0, we have that the inner force is 0
    zeroIdx = find(ts == 0);
    if ~isempty(zeroIdx)
       Fs(zeroIdx) = 0;
    end
    
    % Extract quantities from t > 0
    ds = ds(ts > 0);
    d_ts = d_ts(ts > 0);
    Js = Js(ts > 0);

    % Solve for c(t)
    options = optimoptions('fsolve', 'OptimalityTolerance', 1e-10, ...
        'StepTolerance',1e-10);
    cGuess = 2 * sqrt(pi * ds ./ (4 * epsilon^2 * Js));
    zero_fun = @(c) c.^2 + 4 * c + 2 * log(c) + 1 ...
        - pi * ds ./ (epsilon^2 * Js);
    cs = fsolve(zero_fun, cGuess, options);
    
    % Return forces
    Fs(ts > 0) = (8 * epsilon^4 * d_ts.^2 .* Js.^2 .* cs / pi) ...
        .* (pi * ds ./ (epsilon^2 * Js) - cs.^2 / 3 - 2 * cs - 2 * log(cs) +1);

end