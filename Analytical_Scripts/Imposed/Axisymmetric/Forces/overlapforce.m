function Fs = overlapforce(ds, d_ts, epsilon)
%%innerforce
%   Returns the overlap solution for the force on the substrate in the
%   two-dimensional, plate impact case.
    
    Fs = 16 * sqrt(2) * epsilon * ds.^3 .* d_ts.^2 / 9;
    
    % At t == 0, set the force to be 0
    zeroIdx = find(ds == 0);
    if ~isempty(zeroIdx)
       Fs(zeroIdx) = 0;
    end
end