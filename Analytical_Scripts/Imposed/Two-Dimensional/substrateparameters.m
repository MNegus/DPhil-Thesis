function [epsilon, k, q, omega] = substrateparameters()
%QUADRATICPARAMETERS Set consistent parameters for quadratic substrate

    epsilon = 0.1;
    k = 2 / (3 * epsilon);
    q = 0.1;
	omega = 3.0;
end

