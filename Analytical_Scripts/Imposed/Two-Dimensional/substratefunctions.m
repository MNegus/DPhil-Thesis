function SubstrateFunctions = substratefunctions(type)
%SUBSTRATEFUNCTIONS Returns substrate coefficient functions for the imposed
%substrate.
%   We assume a quadratic substrate of the form
%   w(x, t) = w(t) (1 - x^2 / L^2),
%   where
%   w(t) = q * (t + 1 - cos(omega * t)),
%   we return the substrate coefficients such that
%   w(x, t) = a(t) + b(t) * x^2.
%   type needs to be "stationary" (so a = b = 0), "flat" (so a = w(t) and b
%   = 0) or "curved" (so a = w(t) and b = -k^2 * w(t)).

    %% Check for valid input
    if (type ~= "stationary") && (type ~= "flat") && (type ~= "curved")
        error("Invalid type. Needs to either be 'stationary', 'flat' or 'curved'.");
    end

    %% Load in parameters
    [epsilon, k, q, omega] = substrateparameters();

    %% Set substrate coefficients
    if type == "stationary"
        SubstrateFunctions.a = @(t) zeros(size(t));
        SubstrateFunctions.a_t = @(t) zeros(size(t));
        SubstrateFunctions.a_tt = @(t) zeros(size(t));
    else
        SubstrateFunctions.a = @(t) q * (t.^2 + 1 - cos(omega * t));
        SubstrateFunctions.a_t = @(t) q * (2 * t + omega * sin(omega * t));
        SubstrateFunctions.a_tt = @(t) q * (2 + omega^2 * cos(omega * t));
    end
    
    % Set b coefficients
    if type == "curved"
        SubstrateFunctions.b = @(t) -k^2 * SubstrateFunctions.a(t);
        SubstrateFunctions.b_t = @(t) - k^2 * SubstrateFunctions.a_t(t);
        SubstrateFunctions.b_tt = @(t) - k^2 * SubstrateFunctions.a_tt(t);
    else
        SubstrateFunctions.b = @(t) zeros(size(t));
        SubstrateFunctions.b_t = @(t) zeros(size(t));
        SubstrateFunctions.b_tt = @(t) zeros(size(t));
    end
    
    %% Define outer functions (scaled with epsilons where appropriate)
    % aHat coefficients
    aHat = @(t) SubstrateFunctions.a(t);
    aHat_t = @(t) SubstrateFunctions.a_t(t);
    aHat_tt = @(t) SubstrateFunctions.a_tt(t);
    
    % bHat coefficients
    bHat = @(t) epsilon^2 * SubstrateFunctions.b(t);
    bHat_t = @(t) epsilon^2 * SubstrateFunctions.b_t(t);
    bHat_tt = @(t) epsilon^2 * SubstrateFunctions.b_tt(t);
    
    %% Turnover points
    d = @(t) 2 * sqrt((t - aHat(t)) ./ (1 + 2 * bHat(t)));
    d_t = @(t) ((1 + 2 * bHat(t)) .* (1 - aHat_t(t)) - 2 * (t - aHat(t)) .* bHat_t(t)) ...
        ./ (sqrt(t - aHat(t)) .* (1 + 2 * bHat(t)).^(3/2));
    d_tt = @(t) -((1 + 2 * bHat(t)) .* d_t(t).^2 ...
        + 4 * bHat_t(t) .* d(t) .* d_t(t) ...
        + 2 * aHat_tt(t) ...
        + bHat_tt(t) .* d(t).^2) ./ ((1 + 2 * bHat(t)) .* d(t));
    
    
    %% A, B, C functions
    % B, jet-thickness factor
    B = @(t) aHat_t(t) + 0.5 * bHat_t(t) .* d(t).^2;
    
    % C function
    function Cval = CFullFun(t, d, d_t, B, aHat_t)
       % Cs, overlap pressure factor
        Cval = d(t) .* d_t(t) .* (1 - B(t));
        
        zeroIdx = find(t == 0);
        if ~isempty(zeroIdx)
           Cval(zeroIdx) = 2 * (1 - aHat_t(0)); 
        end 
    end
    C = @(t) CFullFun(t, d, d_t, B, aHat_t);

    % A function
    A = @(t) C(t) - 0.5 * aHat_tt(t) .* d(t).^2 ...
        - (1 / 8) * bHat_tt(t) .* d(t).^4;
    
    %% Jet thickness
    SubstrateFunctions.J = @(t) pi * (1 - B(t)).^2 .* d(t) ./ (8 * d_t(t).^2);
    
    %% Full substrate solution (only works for t being a scalar)
    SubstrateFunctions.w = @(x, t) SubstrateFunctions.a(t) * ones(size(x)) ...
        + SubstrateFunctions.b(t) * x.^2;
    
    %% Load in remaining functions
    SubstrateFunctions.aHat = aHat;
    SubstrateFunctions.aHat_t = aHat_t;
    SubstrateFunctions.aHat_tt = aHat_tt;
    SubstrateFunctions.bHat = bHat;
    SubstrateFunctions.bHat_t = bHat_t;
    SubstrateFunctions.bHat_tt = bHat_tt;
    
    SubstrateFunctions.d = d;
    SubstrateFunctions.d_t = d_t;
    SubstrateFunctions.d_tt = d_tt;
    
    SubstrateFunctions.A = A;
    SubstrateFunctions.B = B;
    SubstrateFunctions.C = C;
    
end

