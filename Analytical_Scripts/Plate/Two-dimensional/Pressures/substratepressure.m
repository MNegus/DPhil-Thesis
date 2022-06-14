function [ps_composite, ps_outer, ps_inner] ...
    = substratepressure(xs, PlatePositions, TimeDependents, epsilon)
%%substratepressure
%   Master function for determining the pressure along the substrate in the
%   two-dimensional, plate impact case. 
%   Inputs:
%   * xs: The x coordinates along the plate
%   * PlatePositions: Structure array for the current plate position
%   variables
%   * TimeDependents: Structure array for the current time-dependent
%   quantities (d, d_t, J, A, B, C).
%   * epsilon: The small time parameters
%
%   Outputs:
%   * ps_composite: The composite solution for the pressure
%   * ps_outer: The outer solution for the pressure
%   * ps_inner: The inner solution for the pressure

    %% Load in quantities from structure arrays
    w_tt = PlatePositions.w_tts;
    d = TimeDependents.ds;
    d_t = TimeDependents.d_ts;
    J = TimeDependents.Js;
    A = TimeDependents.As;
    C = TimeDependents.Cs;
    
    %% Determine pressures
    ps_outer = outerpressure(xs, d, A, w_tt, epsilon);
    ps_inner = innerpressure(xs, d, d_t, J, epsilon);
    ps_overlap = overlappressure(xs, d, C, epsilon);
    ps_composite = compositepressure(ps_outer, ps_inner, ps_overlap);