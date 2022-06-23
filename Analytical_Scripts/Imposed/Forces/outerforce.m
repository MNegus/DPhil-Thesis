function Fs = outerforce(ts, SubstrateFunctions)
%%outerforce
%   Returns the outer solution for the force on the substrate in the
%   two-dimensional, plate impact case.

    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional solution
        Fs = pi * SubstrateFunctions.A(ts);
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric solution
        
        % Load functions
        epsilon = SubstrateFunctions.epsilon;
        ds = SubstrateFunctions.d(ts);
        d_ts = SubstrateFunctions.d_t(ts);
        d_tts = SubstrateFunctions.d_tt(ts);
        
        % Set the force
        Fs = (8 * epsilon / 9) * ds.^3 .* (ds .* d_tts + 4 * d_ts.^2);
    
        % At t == 0, set the force to be 0
        zeroIdx = find(ts == 0);
        if ~isempty(zeroIdx)
           Fs(zeroIdx) = 0;
        end
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end
end