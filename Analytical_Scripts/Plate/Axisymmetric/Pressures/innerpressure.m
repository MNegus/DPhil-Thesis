function ps = innerpressure(rs, d, d_t, J, epsilon)
%%innerpressure
% Returns the inner solution for the pressure along the substrate in the
% two-dimensional, plate impact case. 

    %% Define inner variables
    turnoverIdx = sum(rs < epsilon * d); % Index of the turnover point
    xTildes = (rs - epsilon * d) / epsilon^3;

    %% Determine sigma parameters
    % The xTilde values are found as function of sigma, where
    % xTilde = (J(t)/pi)*(2 sigma - exp(-2 * sigma) - 4 * exp(-sigma)- 1),
    % so we define our own parameter set of sigmas, and then interpolate
    % these back onto xTildes.
    
    % Negative sigmas, where xTilde -> - Infinity
    negative_sigmas = -log(sqrt(4 - pi * xTildes(1 : turnoverIdx) / J) - 2);
    
    % Positive sigmas, where xTilde -> +Infinity
    positive_sigmas = pi * xTildes(turnoverIdx + 1 : end) / (2 * J);
    
    % Combine sigma arrays
    sigmas = cat(2, negative_sigmas, positive_sigmas);
    sigmas = unique(sigmas); % Remove any repeats
    
    % Save xs as a function of sigmas
    xsSigmas = epsilon * d ...
        + (epsilon^3 * J / pi) * (2 * sigmas - 4 * exp(-sigmas) - exp(-2 * sigmas) - 1);
    
    %% Find pressure as a function of sigma
    psSigmas = (1 / epsilon^2) ...
        * 2 * d_t^2 * exp(-sigmas) ./ (1 + exp(-sigmas)).^2;
    
    %% Interpolate the parameterised solution onto the xs dependent solution
    ps = interp1(xsSigmas, psSigmas, rs);
    
end