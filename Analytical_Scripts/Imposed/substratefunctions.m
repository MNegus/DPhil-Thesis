function SubstrateFunctions = substratefunctions(type, dimension)
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
    % Valid substrate type
    if (type ~= "stationary") && (type ~= "flat") && (type ~= "curved")
        error("Invalid type. Needs to either be 'stationary', 'flat' or 'curved'.");
    end
    
    % Valid dimension
    if (dimension ~= "2D") && (dimension ~= "axi")
        error("Invalid dimension. Needs to either be '2D' or 'axi'.");
    end
    
    % Error if try curved and axisymmetric
    if (dimension == "axi") && (type == "curved")
        error("Cannot have a curved substrate in axisymmetric dimension.");
    end

    %% Load in parameters
    [epsilon, k, q, omega] = substrateparameters();

    %% Set parameters
    SubstrateFunctions.type = type;
    SubstrateFunctions.dimension = dimension;
    SubstrateFunctions.epsilon = epsilon;
    SubstrateFunctions.k = k;
    SubstrateFunctions.q = q;
    SubstrateFunctions.omega = omega;
    
    %% Set substrate coefficients
    if type == "stationary"
        a = @(t) zeros(size(t));
        
        a_t = @(t) zeros(size(t));
        
        a_tt = @(t) zeros(size(t));
    else
        a = @(t) q * (t.^2 + 1 - cos(omega * t));
        
        a_t = @(t) q * (2 * t + omega * sin(omega * t));
        
        a_tt = @(t) q * (2 + omega^2 * cos(omega * t));
    end
    
    % Set b coefficients
    if type == "curved"
%         b = @(t) -k^2 * a(t);
%         b_t = @(t) - k^2 * a_t(t);
%         b_tt = @(t) - k^2 * a_tt(t);

        % Alt: Curve crosses x axis at c = epsilon * p * sqrt(t) 
        p = 1.25;
        b = @(t) -a(t) ./ (p^2 * t);
        b_t = @(t) (1 / (p^2)) * (a(t) ./ t.^2 - a_t(t) ./ t);
        b_tt = @(t) (1 / (p^2)) * (-2 * a(t) ./ t.^3 + 2 * a_t(t) ./ t.^2 ...
              - a_tt(t) ./ t);
    else
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
        %% Turnover points
        d = @(t) 2 * sqrt((t - a(t)) ./ (1 + 2 * b(t)));
        
        d_t = @(t) ((1 + 2 * b(t)) .* (1 - a_t(t)) - 2 * (t - a(t)) .* b_t(t)) ...
            ./ (sqrt(t - a(t)) .* (1 + 2 * b(t)).^(3/2));
        d_tt = @(t) -((1 + 2 * b(t)) .* d_t(t).^2 ...
            + 4 * b_t(t) .* d(t) .* d_t(t) ...
            + 2 * a_tt(t) ...
            + b_tt(t) .* d(t).^2) ./ ((1 + 2 * b(t)) .* d(t));

        %% A, B, C functions
        % B, jet-thickness factor
        B = @(t) BFullFun(t, d, a_t, b_t);
        SubstrateFunctions.B = B;
    
        C = @(t) CFullFun(t, d, d_t, B, a_t);
        SubstrateFunctions.C = C;
        
        % A function
        A = @(t) C(t) - 0.5 * a_tt(t) .* d(t).^2 ...
            - (1 / 8) * b_tt(t) .* d(t).^4;
        SubstrateFunctions.A = A;

        %% Jet thickness
        SubstrateFunctions.J = @(t) pi * (1 - B(t)).^2 .* d(t) ./ (8 * d_t(t).^2);

        %% Full substrate solution (only works for t being a scalar)
        SubstrateFunctions.w = @(x, t) a(t) * ones(size(x)) ...
            + b(t) * x.^2 / epsilon^2;
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
        
        %% Turnover curves
        d = @(t) sqrt(3 * (t - w(t)));
        d_t = @(t) sqrt(3) * (1 - w_t(t)) ./ (2 * sqrt(t - w(t)));
        d_tt = @(t) - sqrt(3) * ((1 - w_t(t)).^2 + 2 * (t - w(t)) .* w_tt(t)) ...
            ./ (4 * (t - w(t)).^(3/2));
        
        %% Jet thickness
        SubstrateFunctions.J = @(t) 2 * d(t).^3 / (9 * pi);
        
        %% Pressure coefficients
        SubstrateFunctions.A = @(t) (4 * d(t).^3 / 9) .* (4 * d_t(t).^2 + d(t) .* d_tt(t));
        SubstrateFunctions.C = @(t) 2 * sqrt(2) * d(t).^(3/2) * d_t(t).^2 / (3 * pi);
    end
    
    %% Load in turnover point/curve functions
    SubstrateFunctions.d = d;
    SubstrateFunctions.d_t = d_t;
    SubstrateFunctions.d_tt = d_tt;
    
    %% Function definitions
    % B function
    function Bval = BFullFun(t, d, aHat_t, bHat_t)
        Bval = aHat_t(t) + 0.5 * bHat_t(t) .* d(t).^2;
        
        zeroIdx = find(t == 0);
        if ~isempty(zeroIdx)
           Bval(zeroIdx) = aHat_t(0); 
        end 
    end

    % C function
    function Cval = CFullFun(t, d, d_t, B, aHat_t)
       % Cs, overlap pressure factor
        Cval = d(t) .* d_t(t) .* (1 - B(t));

        zeroIdx = find(t == 0);
        if ~isempty(zeroIdx)
           Cval(zeroIdx) = 2 * (1 - aHat_t(0)); 
        end 
    end
end

