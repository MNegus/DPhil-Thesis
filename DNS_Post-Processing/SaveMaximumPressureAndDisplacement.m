%% SaveMaximumPressure

% Parent directory
% parent_dir = "/home/michael/scratch/DPhil_DNS_Data/Stationary_Membrane";
parent_dir ="/home/michael/scratch/DPhil_DNS_Data/RubberRuns/ALPHA_1.1_GAMMA_668";

% Number of outputs
noOutputs = 7001;

% Maximum pressure and displacement values
wMaxs = zeros(1, noOutputs);
w_tMaxs = zeros(1, noOutputs);
pMaxs = zeros(1, noOutputs);
p0s = zeros(1, noOutputs);

% Loop and save outputs
for k = 1 : noOutputs
    k

    % Load membrane array
    membrane_mat ...
        = importdata(sprintf("%s/membrane_outputs/membrane_arr_%d.txt", ...
            parent_dir, k - 1));

    % Save maximum values
    wMaxs(k) = max(membrane_mat(:, 2));
    w_tMaxs(k) = max(membrane_mat(:, 3));
    pMaxs(k) = max(membrane_mat(:, 4));

    % Save pressure at origin
    x0Idx = find(membrane_mat(:, 1) == 0);
    p0s(k) = membrane_mat(x0Idx, 4);
end

%% Save outputs
MaxStruct.wMaxs = wMaxs;
MaxStruct.w_tMaxs = w_tMaxs;
MaxStruct.pMaxs = pMaxs;
MaxStruct.p0s = p0s;

save(sprintf("%s/MaxStruct.mat", parent_dir), "MaxStruct");
