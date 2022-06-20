function [Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(ts, PlatePositions, TimeDependents, epsilon)
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

    %% Load in quantities from structure arrays
    w_tts = PlatePositions.w_tts;
    ds = TimeDependents.ds;
    d_ts = TimeDependents.d_ts;
    Js = TimeDependents.Js;
    As = TimeDependents.As;
    Cs = TimeDependents.Cs;
    
    %% Determines pressures
    Fs_outer = outerforce(ds, As, w_tts);
    Fs_inner = innerforce(ts, ds, d_ts, Js, Cs, epsilon);
    Fs_overlap = overlapforce(Cs);
    Fs_composite = compositeforce(Fs_outer, Fs_inner, Fs_overlap);

end