function [epsilon, q, omega, p, L] = substrateparameters(type)
%QUADRATICPARAMETERS Set consistent parameters for moving substrate.
%   The substrate parameters are defined such that we have for the flat
%   substrate,
%   w(x, t) = w(t) = q * (t.^2 + 1 - cos(omega * t)),
%   and for the curved substrate,
%   w(x, t) = w(t) (1 - k^2 * x.^2).
%   
%   The type parameter denotes the type of imposed substrate we have. If
%   "type" is "stationary", "flat" or "curved", then we have the parameters
%   from the analytical chapter. If "type" is "flatDNS" or "curvedDNS",
%   then we have the parameters for the DNS chapter.

    if ~exist('type','var')
        % If type is not given, default it to "stationary"
        type = "stationary";
    end

    if (type == "stationary") || (type == "flat") || (type == "curved") 
        %% Analytical chapter parameters
        % Small time parameter
        epsilon = 0.1;

        % Oscillation magnitude
        q = 0.1;

        % Oscillation frequency
        omega = 3.0;

        % Curved substrate crossed x axis at c = epsilon * p * sqrt(t)
        p = 1.25; 

        % Substrate radius
        L = 2;
    elseif (type == "flatDNS")
        %% Flat DNS validation parameters
        % Small time parameter (O(1) in this case)
        epsilon = 1;
        
        % Oscillation magnitude
        q = 0.00125;
        
        % Oscillation frequency
        omega = 12;
        
        % Curved susbtrate parameter (not used here)
        p = 1;
        
        % Substrate radius
        L = 2;
    elseif (type == "curvedDNS")
        %% Curved DNS validation parameters
        % Small time parameter (O(1) in this case)
        epsilon = 1;
        
        % Oscillation magnitude
        q = 0.05;
        
        % Oscillation frequency
        omega = 4;
        
        % Curved susbtrate parameter (not used here)
        p = 1;
        
        % Substrate radius
        L = 2;
    else
        error("Invalid type.");
    end
end

