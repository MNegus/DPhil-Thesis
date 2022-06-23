function [ps_composite, ps_outer, ps_inner] ...
    = substratepressure(rs, TimeDependents, epsilon)
%%substratepressure
%   Master function for determining the pressure along the substrate in the
%   axisymmetric, plate case. 
%   Inputs:
%   * rs: The r coordinates along the plate
%   * TimeDependents: Structure array for the current time-dependent
%   quantities (d, d_t, d_tt, J).
%   * epsilon: The small time parameters
%
%   Outputs:
%   * ps_composite: The composite solution for the pressure
%   * ps_outer: The outer solution for the pressure
%   * ps_inner: The inner solution for the pressure

    %% Load in quantities from structure arrays
    d = TimeDependents.ds;
    d_t = TimeDependents.d_ts;
    d_tt = TimeDependents.d_tts;
    J = TimeDependents.Js;
    
    %% Determine pressures
    ps_outer = outerpressure(rs, d, d_t, d_tt, epsilon);
    ps_inner = innerpressure(rs, d, d_t, J, epsilon);
    ps_overlap = overlappressure(rs, d, d_t, epsilon);
    ps_composite = compositepressure(ps_outer, ps_inner, ps_overlap);
end