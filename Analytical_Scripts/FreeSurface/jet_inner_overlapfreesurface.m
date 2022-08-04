function hs = jet_inner_overlapfreesurface(xs, t, SubstrateFunctions)
%OVERLAPFREESURFACE Summary of this function goes here
%   Detailed explanation goes here

    if SubstrateFunctions.dimension == "2D"
        error("Overlap free surface only supported for dimension == 'axi' currently");
    else
        epsilon = SubstrateFunctions.epsilon;
        w = SubstrateFunctions.w(t);
        d = SubstrateFunctions.d(t);
        J = SubstrateFunctions.J(t);

        xTildes = (xs - epsilon * d) / epsilon^3;
        hs = - epsilon^2 * w ...
            + J * epsilon^3 * (1 + (4 / pi) * exp(-pi * xTildes / (2 * J)));
    end
end
