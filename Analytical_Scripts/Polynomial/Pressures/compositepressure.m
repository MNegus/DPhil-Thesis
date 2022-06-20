function ps = compositepressure(ps_outer, ps_inner, ps_overlap)
%%compositepressure
%   Returns the additive composite expansion of the pressure along the
%   substrate.
    ps = ps_outer + ps_inner - ps_overlap;
end