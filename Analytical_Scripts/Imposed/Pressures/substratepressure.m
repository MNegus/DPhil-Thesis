function [ps_composite, ps_outer, ps_inner] ...
    = substratepressure(xs, t, SubstrateFunctions)
%%substratepressure
%   Master function for determining the pressure along the substrate in the
%   two-dimensional, quadratic substrate case. 
%   Inputs:
%   * xs: The x coordinates along the plate
%   * SubstrateCoefficients: Structure array for the coefficients of the
%   substrate, such that w(x, t) = a + b * x^2
%   * TimeDependents: Structure array for the current time-dependent
%   quantities (d, d_t, J, A, B, C).
%   * epsilon: The small time parameters
%
%   Outputs:
%   * ps_composite: The composite solution for the pressure
%   * ps_outer: The outer solution for the pressure
%   * ps_inner: The inner solution for the pressure
    
    %% Determine pressures
    ps_outer = outerpressure(xs, t, SubstrateFunctions);
    ps_inner = innerpressure(xs, t, SubstrateFunctions);
    ps_overlap = overlappressure(xs, t, SubstrateFunctions);
    ps_composite = compositepressure(ps_outer, ps_inner, ps_overlap);
end