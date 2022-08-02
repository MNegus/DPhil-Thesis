function SubstrateFunctions = imposedsubstratefunctions(type, dimension)
%IMPOSEDSUBSTRATEFUNCTIONS Returns substrate coefficient functions for the 
% imposed substrate. The SubstrateFunctions struct is then fed into
% substratedependents which will solve for the substrate-dependent
% quantities such as the turnover point and jet thickness.

    %% Check for valid input
    % Valid substrate type
    if (type ~= "stationary") && (type ~= "flat") && (type ~= "curved") ...
            && (type ~= "flatDNS") && (type ~= "curvedDNS")
        error("Invalid type. Needs to either be 'stationary', 'flat', 'curved', 'flatDNS' or 'curvedDNS'.");
    end
    
    % Valid dimension
    if (dimension ~= "2D") && (dimension ~= "axi")
        error("Invalid dimension. Needs to either be '2D' or 'axi'.");
    end
    
    % Error if try curved and axisymmetric
    if (dimension == "axi") && ((type == "curved") || (type == "curvedDNS"))
        error("Cannot have a curved substrate in axisymmetric dimension.");
    end

    %% Load in parameters
    [epsilon, q, omega, p, L] = substrateparameters(type);

    %% Set parameters
    SubstrateFunctions.type = type; % Substrate type
    SubstrateFunctions.dimension = dimension; % Dimension
    SubstrateFunctions.epsilon = epsilon; % Small time parameter
    SubstrateFunctions.q = q; % Magnitude of imposed substrate
    SubstrateFunctions.omega = omega; % Angular velocity of imposed substrate
    SubstrateFunctions.p = p; % Such that substrate crosses x axis at epsilon * p * sqrt(t) 
    SubstrateFunctions.L = L; % Substrate radius
    
    %% Set substrate coefficients
    if type == "stationary"
        a = @(t) zeros(size(t));
        
        a_t = @(t) zeros(size(t));
        
        a_tt = @(t) zeros(size(t));
    else
        % Imposed susbtrate motion
        a = @(t) q * (t.^2 + 1 - cos(omega * t));
        
        a_t = @(t) q * (2 * t + omega * sin(omega * t));
        
        a_tt = @(t) q * (2 + omega^2 * cos(omega * t));
    end
    
    % Set b coefficients
    if type == "curved"
        % Curved substrate crosses x axis at c = epsilon * p * sqrt(t) 
        p = 1.25;
        b = @(t) -a(t) ./ (p^2 * t);
        b_t = @(t) (1 / (p^2)) * (a(t) ./ t.^2 - a_t(t) ./ t);
        b_tt = @(t) (1 / (p^2)) * (-2 * a(t) ./ t.^3 + 2 * a_t(t) ./ t.^2 ...
              - a_tt(t) ./ t);
    elseif type == "curvedDNS"
        % Curved DNS substrate crosses x axis at x = L
        b = @(t) -epsilon^2 * a(t) / L^2;
        b_t = @(t) -epsilon^2 * a_t(t) / L^2;
        b_tt = @(t) -epsilon^2 * a_tt(t) / L^2;
    else
        % Else we have no curved part
        b = @(t) zeros(size(t));
        b_t = @(t) zeros(size(t));
        b_tt = @(t) zeros(size(t));
    end
    
    % Load functions into structure array
    SubstrateFunctions.a = a;
    SubstrateFunctions.a_t = a_t;
    SubstrateFunctions.a_tt = a_tt;
    
    SubstrateFunctions.b = b;
    SubstrateFunctions.b_t = b_t;
    SubstrateFunctions.b_tt = b_tt;
    
    if dimension == "2D"
    %% Two-dimensional quantities
        %% Full substrate solution (only works for t being a scalar)
        SubstrateFunctions.w = @(x, t) a(t) * ones(size(x)) ...
            + b(t) * x.^2 / epsilon^2;
        SubstrateFunctions.w_t = @(x, t) a_t(t) * ones(size(x)) ...
            + b_t(t) .* x.^2 / epsilon^2;
        SubstrateFunctions.w_xx = @(x, t) 2 * b(t) / epsilon^2;
        SubstrateFunctions.w_xt = @(x, t) 2 * b_t(t) * x / epsilon^2;
        SubstrateFunctions.w_tt = @(x, t) a_tt(t) + b_tt(t) * x.^2 / epsilon^2;
        
    else
    %% Axisymmetric quantities
        %% Substrate positions
        w = SubstrateFunctions.a;
        SubstrateFunctions.w = w;
        
        w_t = SubstrateFunctions.a_t;
        SubstrateFunctions.w_t = w_t;
        
        w_tt = SubstrateFunctions.a_tt;
        SubstrateFunctions.w_tt = w_tt;
    end
    
    
    %% Load in substrate dependent functions
    SubstrateFunctions = substratedependents(SubstrateFunctions);
    
end

