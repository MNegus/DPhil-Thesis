function newSubstrateFunctions = substratedependents(SubstrateFunctions)
%SUBSTRATEDEPENDENTS Determines substrate-dependent functions
%   The input SubstrateFunctions should have the dimension (2D or axi), and
%   the substrate functions stored in it, i.e. functions a(t) and b(t) such
%   that w(x, t) = a(t) + b(t) * x.^2, as well as the derivatives of a and
%   b. Then, substrate dependent quantities (such as the turnover point
%   location and jet thickness) are added to the struct. 


    %% Define new struct to fill
    newSubstrateFunctions = SubstrateFunctions;

    %% Load in useful variables and functions
    dimension = SubstrateFunctions.dimension; % Dimension (either 2D or axi)
    
    % Rigid term
    a = SubstrateFunctions.a;
    a_t = SubstrateFunctions.a_t;
    a_tt = SubstrateFunctions.a_tt;
    
    % Curved term
    b = SubstrateFunctions.b;
    b_t = SubstrateFunctions.b_t;
    b_tt = SubstrateFunctions.b_tt;


    %% Solutions depending on dimension
    if dimension == "2D"
    %% Two-dimensional quantities (for quadratic substrate)
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
        newSubstrateFunctions.B = B;
    
        C = @(t) CFullFun(t, d, d_t, B, a_t);
        newSubstrateFunctions.C = C;
        
        % A function
        A = @(t) C(t) - 0.5 * a_tt(t) .* d(t).^2 ...
            - (1 / 8) * b_tt(t) .* d(t).^4;
        newSubstrateFunctions.A = A;

        %% Jet thickness
        newSubstrateFunctions.J = @(t) pi * (1 - B(t)).^2 .* d(t) ./ (8 * d_t(t).^2);

    else
    %% Axisymmetric quantities (for flat substrate)
        
        %% Turnover curves
        d = @(t) sqrt(3 * (t - a(t)));
        d_t = @(t) sqrt(3) * (1 - a_t(t)) ./ (2 * sqrt(t - a(t)));
        d_tt = @(t) - sqrt(3) * ((1 - a_t(t)).^2 + 2 * (t - a(t)) .* a_tt(t)) ...
            ./ (4 * (t - a(t)).^(3/2));
        
        %% Jet thickness
        newSubstrateFunctions.J = @(t) 2 * d(t).^3 / (9 * pi);
        
        %% Pressure coefficients
        newSubstrateFunctions.A = @(t) (4 * d(t).^3 / 9) .* (4 * d_t(t).^2 + d(t) .* d_tt(t));
        newSubstrateFunctions.C = @(t) 2 * sqrt(2) * d(t).^(3/2) * d_t(t).^2 / (3 * pi);
    end
    
    %% Load in turnover point/curve functions
    newSubstrateFunctions.d = d;
    newSubstrateFunctions.d_t = d_t;
    newSubstrateFunctions.d_tt = d_tt;
    
    %% Function definitions
    % B function
    function Bval = BFullFun(t, d, a_t, b_t)
        Bval = a_t(t) + 0.5 * b_t(t) .* d(t).^2;
        
        zeroIdx = find(t == 0);
        if ~isempty(zeroIdx)
           Bval(zeroIdx) = a_t(0); 
        end 
    end

    % C function
    function Cval = CFullFun(t, d, d_t, B, a_t)
       % Cs, overlap pressure factor
        Cval = d(t) .* d_t(t) .* (1 - B(t));

        zeroIdx = find(t == 0);
        if ~isempty(zeroIdx)
           Cval(zeroIdx) = 2 * (1 - a_t(0)); 
        end 
    end
end