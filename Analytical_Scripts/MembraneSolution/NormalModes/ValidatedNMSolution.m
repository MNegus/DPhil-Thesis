function [N, delta_d, ts, NM_Struct] ...
    = ValidatedNMSolution(alpha, beta, gamma, epsilon, L, tmax, delta_t)
%VALIDATEDNMSOLUTION Finds the normal modes solution with the appropriate N

    %% Derived parameters
    ts = 0 : delta_t : tmax;

    %% Parameters to be passed in
    q = 10;
    tol = 1e-4;
    N_MEMBRANE = 21848;
    DELTA_X = L / (N_MEMBRANE - 1); 
    xs = (0 : DELTA_X : L - DELTA_X)';

    %% d dependent parameters
    d_max = 2 * sqrt(tmax);

    %% Initial guess for delta_d and N
    delta_d = delta_t / sqrt(tmax);

    %% Loops until found appropriate delta_d and N
    converged = 0;
    while (converged == 0)
        %% Determines N_max for current delta_d
        delta_d
        N_max = NMax(alpha, beta, gamma, L, q, epsilon, delta_t)
        
        %% Initialise N and solves ode
        N0 = min(64, floor(N_max / 4));
        N = N0;
        NM_Struct_d_form ...
            = NormalModesODE(alpha, beta, gamma, epsilon, delta_d, d_max, N, L);

        % Converts solution to t-form
        NM_Struct = NormalModesTemporalForm(ts, NM_Struct_d_form, ...
            alpha, epsilon, delta_t);
        
        % Determines tolerance depending on the maximum displacement
        [ws, w_ts, ~] = MembraneSolutionNM(xs, NM_Struct.as(end, :), ...
            NM_Struct.a_ts(end, :), NM_Struct.q_ts(end, :), ...
            NM_Struct.ds(end), L, N, epsilon);
        tol = min(tol, 0.01 * max(ws))
        
        %% Loops until we've found an N or N >= N_max
        diff = 1e6;
        while ((diff > tol) && (2 * N < N_max))
            % Doubles N
            new_N = 2 * N
            pause(1)

            % Solves ode with new N
            new_NM_Struct_d_form ...
                = NormalModesODE(alpha, beta, gamma, epsilon, delta_d, d_max, new_N, L);
            
            % New solution in t form
            new_NM_Struct ...
                = NormalModesTemporalForm(ts, new_NM_Struct_d_form, ...
                    alpha, epsilon, delta_t);
            
            % Compares ws solution for all time between the two
            diff = 0;
%             figure(1);
            for k = length(ts)
                [ws, w_ts, ~] ...
                    = MembraneSolutionNM(xs, NM_Struct.as(k, :), ...
                    NM_Struct.a_ts(k, :), NM_Struct.q_ts(k, :), ...
                    NM_Struct.ds(k), L, N, epsilon);
                [new_ws, new_w_ts, ~] ...
                    = MembraneSolutionNM(xs, new_NM_Struct.as(k, :), ...
                    new_NM_Struct.a_ts(k, :), new_NM_Struct.q_ts(k, :), ...
                    new_NM_Struct.ds(k), L, new_N, epsilon);
                
                norm = max(abs(ws - new_ws));
%                 norm = max(abs((ws - new_ws) ./ new_ws))

                diff = max(diff, norm);
%                 diff = max(diff, max(abs(w_ts - new_w_ts)))
                
                plot(xs, ws);
                hold on;
                plot(xs, new_ws);
                hold off;
                title(['k = ', num2str(k), ', new_N = ', num2str(new_N), ', delta_d = ', num2str(delta_d)]);
                drawnow;
                pause(0.0000001);

                if (diff > tol)
                    delta_d
                    N
                    break;
                end
                
            end
            
            % Update solutions
            N = new_N;
            NM_Struct = new_NM_Struct;
        end
        
        %% Checks for value of N
        if (diff < tol)
            converged = 1;
        else
%            N = N0;
           delta_t = delta_t / 10;
           delta_d = delta_t / sqrt(tmax);
        end
    
        
    end

end