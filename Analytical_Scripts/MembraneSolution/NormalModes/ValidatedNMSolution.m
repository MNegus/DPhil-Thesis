function [N, delta_d, ts, ds, as, a_ts, a_tts, q_ts, d_ts, Js] ...
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
        [t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, d_ts_d_form, Js_d_form, omegas] ...
            = NormalModesODE(alpha, beta, gamma, epsilon, delta_d, d_max, N, L);
        % Converts solution to t-form
        [ds, as, a_ts, a_tts, q_ts, d_ts, Js] = NormalModesTemporalForm(ts, ...
            t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, d_ts_d_form, Js_d_form, ...
            omegas, alpha, epsilon, delta_t);
        
        % Determines tolerance depending on the maximum displacement
        [ws, w_ts, ~] = MembraneSolutionNM(xs, as(end, :), a_ts(end, :), q_ts(end, :), ds(end), L, N, epsilon);
        tol = min(tol, 0.1 * max(ws))
        
        %% Loops until we've found an N or N >= N_max
        diff = 1e6;
        while ((diff > tol) && (2 * N < N_max))
            % Doubles N
            new_N = 2 * N
            pause(1)

            % Solves ode with new N
            [new_t_vals_d_form, new_d_vals_d_form, new_as_d_form, ...
                new_a_ts_d_form, new_d_ts_d_form, new_Js_d_form, new_omegas] ...
                = NormalModesODE(alpha, beta, gamma, epsilon, delta_d, d_max, new_N, L);
            
            % New solution in t form
            [new_ds, new_as, new_a_ts, new_a_tts, new_q_ts, new_d_ts, new_Js] ...
                = NormalModesTemporalForm(ts, new_t_vals_d_form, ...
                    new_d_vals_d_form, new_as_d_form, new_a_ts_d_form, ...
                    new_d_ts_d_form, new_Js_d_form, new_omegas, ...
                    alpha, epsilon, delta_t);
            
            % Compares ws solution for all time between the two
            diff = 0;
%             figure(1);
            for k = length(ts)
                [ws, w_ts, ~] = MembraneSolutionNM(xs, as(k, :), a_ts(k, :), q_ts(k, :), ds(k), L, N, epsilon);
                [new_ws, new_w_ts, ~] = MembraneSolutionNM(xs, new_as(k, :), new_a_ts(k, :), new_q_ts(k, :), new_ds(k), L, new_N, epsilon);
                
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
            ds = new_ds;
            as = new_as;
            a_ts = new_a_ts;
            a_tts = new_a_tts;
            q_ts = new_q_ts;
            d_ts = new_d_ts;
            Js = new_Js;
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