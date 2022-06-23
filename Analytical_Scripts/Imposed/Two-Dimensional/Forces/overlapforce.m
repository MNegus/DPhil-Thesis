function Fs = overlapforce(ts, SubstrateFunctions)
%%innerforce
%   Returns the overlap solution for the force on the substrate in the
%   two-dimensional, plate impact case.
    
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional solution
        Fs = 2 * sqrt(2) * SubstrateFunctions.C(ts);
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric solution
        
        % Load functions
        epsilon = SubstrateFunctions.epsilon;
        ds = SubstrateFunctions.d(ts);
        d_ts = SubstrateFunctions.d_t(ts);
        
        % Set the force
        Fs = 16 * sqrt(2) * epsilon * ds.^3 .* d_ts.^2 / 9;
    
        % At t == 0, set the force to be 0
        zeroIdx = find(ds == 0);
        if ~isempty(zeroIdx)
           Fs(zeroIdx) = 0;
        end
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end
end