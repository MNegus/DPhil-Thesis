function [Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(ts, SubstrateFunctions)
%%substrateforce
%   Master function for determining the force along the substrate in the
%   two-dimensional, plate impact case. 
%   Inputs:
%   * PlatePositions: Structure array for the plate position
%   variables
%   * TimeDependents: Structure array for the time-dependent
%   quantities (d, d_t, J, A, B, C).
%   * epsilon: The small time parameters
%
%   Outputs:
%   * Fs_composite: The composite solution for the force
%   * Fs_outer: The outer solution for the force
%   * Fs_inner: The inner solution for the force

    
    %% Determines forces
    Fs_outer = outerforce(ts, SubstrateFunctions);
    Fs_inner = innerforce(ts, SubstrateFunctions);
    Fs_overlap = overlapforce(ts, SubstrateFunctions);
    Fs_composite = compositeforce(Fs_outer, Fs_inner, Fs_overlap);

end