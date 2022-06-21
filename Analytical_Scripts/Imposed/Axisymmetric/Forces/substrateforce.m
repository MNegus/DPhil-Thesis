function [Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(ts, TimeDependents, epsilon)
%%substrateforce
%   Master function for determining the force along the substrate in the
%   axisymmetric, plate impact case. 
%   Inputs:
%   * PlatePositions: Structure array for the plate position
%   variables
%   * TimeDependents: Structure array for the time-dependent
%   quantities (d, d_t, d_tt, J).
%   * epsilon: The small time parameters
%
%   Outputs:
%   * Fs_composite: The composite solution for the force
%   * Fs_outer: The outer solution for the force
%   * Fs_inner: The inner solution for the force

    %% Load in quantities from structure arrays
    ds = TimeDependents.ds;
    d_ts = TimeDependents.d_ts;
    d_tts = TimeDependents.d_tts;
    Js = TimeDependents.Js;
    
    %% Determines pressures
    Fs_outer = outerforce(ds, d_ts, d_tts, epsilon);
    Fs_inner = innerforce(ts, ds, d_ts, Js, epsilon);
    Fs_overlap = overlapforce(ds, d_ts, epsilon);
    Fs_composite = compositeforce(Fs_outer, Fs_inner, Fs_overlap);

end