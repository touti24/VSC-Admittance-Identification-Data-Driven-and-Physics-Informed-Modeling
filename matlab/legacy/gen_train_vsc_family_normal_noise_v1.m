%% generate_training_data_vsc_family_normal_noise.m
% =========================================================================
% VSC Parameter Family Training Dataset Generator (Normal Noise Regime)
%
% Description:
%   Generates a synthetic training dataset comprising 3,000 randomized
%   small-signal 6-state Voltage Source Converter (VSC) systems.
%
% Methodology:
%   - Parameter Sampling: Physical parameters (filter inductance Lf,
%     filter resistance Rf, current-loop PI gains, PLL PI gains) are sampled
%     within +/-50% of the nominal VSC I design point across active/reactive
%     current operating ranges.
%   - Perturbation Protocol: Two-event orthogonal grid voltage steps
%     (Vd step and Vq step) in the global synchronous dq-frame.
%   - Sensor Noise Modeling: Gaussian measurement noise is embedded at
%     realistic industrial sensor noise levels:
%       * 30 dB (~3.0% relative noise)
%       * 28 dB (~4.0% relative noise)
%       * 26 dB (~5.0% relative noise)
%       * Inf dB (noise-free analytical baseline)
%   - Ground-Truth Calculation: Analytical 2x2 dq-admittance matrices Y(s)
%     are computed on a dense resonance-focused frequency grid.
%
% Output:
%   - Noisy & clean time-domain step-response CSVs per system
%   - Ground-truth admittance matrices and physical parameters (.mat)
%   - Manifest file ('manifest.csv') mapping all systems
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & SIMULATION PARAMETERS
n_systems     = 3000;                       % Target number of stable systems
Fs            = 2500;                       % [Hz] Sampling frequency
T_window      = 0.98;                       % [s] Simulation duration per step event
step_size_V   = 1.0;                        % [V] Small-signal voltage step perturbation
SNR_dB_list   = [30, 28, 26];               % Target SNR levels [dB] (~3% to 5% noise)
out_dir       = "training_data_out_normal_noise";
rng(7);                                     % Enforce exact reproducibility

% VSC I nominal baseline parameters
Lf_nom        = 1.0e-3;                     % [H] Nominal filter inductance
Rf_nom        = 0.3e-3;                     % [Ohm] Nominal filter resistance
Kp_i_nom      = 10.5;                       % Current loop proportional gain
Ki_i_nom      = 5900.0;                     % Current loop integral gain
Kp_pll_nom    = 1.2;                        % SRF-PLL proportional gain
Ki_pll_nom    = 257.0;                      % SRF-PLL integral gain

range_frac    = 0.5;                        % Fractional parameter spread (+/-50%)
Id_range      = [0, 100];                   % [A] Active current operating range
Iq_range      = [-50, 50];                  % [A] Reactive current operating range
Vd_ss         = 311.0;                      % [V] Steady-state peak phase voltage
w0            = 2 * pi * 50;                % [rad/s] Fundamental grid frequency (50 Hz)

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Time vector (column vector)
dt = 1.0 / Fs;
t  = (0 : dt : T_window - dt)';
Nt = length(t);

% Dense non-uniform frequency grid focused on the sub-synchronous resonance band
f_coarse_lo = linspace(1, 35, 60);
f_dense     = linspace(35, 75, 2000);
f_coarse_hi = linspace(75, 100, 60);
freqs_Hz    = unique([f_coarse_lo, f_dense, f_coarse_hi]);

manifest = {};
sys_id   = 0;
attempts = 0;
max_attempts = n_systems * 3;

%% 2. DATASET GENERATION LOOP
fprintf("=========================================================================\n");
fprintf("Starting VSC Training Data Generation (%d Systems, SNR in [30, 28, 26] dB)\n", n_systems);
fprintf("=========================================================================\n");

while sys_id < n_systems && attempts < max_attempts
    attempts = attempts + 1;

    % Random parameter draw around nominal converter family
    Lf     = Lf_nom     * (1 + range_frac * (2 * rand() - 1));
    Rf     = Rf_nom     * (1 + range_frac * (2 * rand() - 1));
    Kp_i   = Kp_i_nom   * (1 + range_frac * (2 * rand() - 1));
    Ki_i   = Ki_i_nom   * (1 + range_frac * (2 * rand() - 1));
    Kp_pll = Kp_pll_nom * (1 + range_frac * (2 * rand() - 1));
    Ki_pll = Ki_pll_nom * (1 + range_frac * (2 * rand() - 1));
    Id_eq  = Id_range(1) + diff(Id_range) * rand();
    Iq_eq  = Iq_range(1) + diff(Iq_range) * rand();

    % Linearized 6-state VSC small-signal model
    % States : [id, iq, xi_d, xi_q, theta_pll, xi_pll]
    % Inputs : [Vd, Vq]
    % Outputs: [id, iq]
    n = 6;
    A = zeros(n, n);
    B = zeros(n, 2);
    C = zeros(2, n);
    D = zeros(2, 2);

    A(1, 1) = -Rf / Lf;      A(1, 2) =  w0;
    A(2, 1) = -w0;           A(2, 2) = -Rf / Lf;
    A(1, 3) =  Kp_i / Lf;    A(2, 4) =  Kp_i / Lf;
    A(3, 1) = -1.0;          A(4, 2) = -1.0;

    A(5, 5) = -Kp_pll * Vd_ss;
    A(5, 6) =  1.0;
    A(6, 5) = -Ki_pll * Vd_ss;
    A(1, 5) =  Iq_eq * w0;
    A(2, 5) = -Id_eq * w0;

    B(1, 1) = -1.0 / Lf;
    B(2, 2) = -1.0 / Lf;
    B(5, 2) =  1.0;

    C(1, 1) = 1.0;
    C(2, 2) = 1.0;

    % Reject unstable parameter combinations
    if any(real(eig(A)) >= 0)
        continue;
    end
    assert(isreal(A) && isreal(B) && isreal(C), "State-space matrices must be real-valued.");
    sysL = ss(A, B, C, D);

    sys_id = sys_id + 1;

    % Two-event voltage step responses
    u_ev1 = zeros(Nt, 2); u_ev1(t >= 0.1, 1) = step_size_V;
    u_ev2 = zeros(Nt, 2); u_ev2(t >= 0.1, 2) = step_size_V;
    y_ev1_clean = lsim(sysL, u_ev1, t);
    y_ev2_clean = lsim(sysL, u_ev2, t);

    % Frequency-domain ground truth Y(s) = C*(sI - A)^(-1)*B + D
    Y_true = zeros(2, 2, length(freqs_Hz));
    for k = 1:length(freqs_Hz)
        s = 1j * 2 * pi * freqs_Hz(k);
        Y_true(:, :, k) = C * ((s * eye(n) - A) \ B) + D;
    end

    % Export clean reference CSV (SNR = Inf)
    T_clean = table(t, y_ev1_clean(:, 1), y_ev1_clean(:, 2), y_ev2_clean(:, 1), y_ev2_clean(:, 2), ...
        'VariableNames', {'time', 'id_ev1', 'iq_ev1', 'id_ev2', 'iq_ev2'});
    writetable(T_clean, fullfile(out_dir, sprintf('sys_%05d_SNR_Inf.csv', sys_id)));

    % Export noisy CSVs across target SNR levels
    for snr_idx = 1:length(SNR_dB_list)
        SNR_dB = SNR_dB_list(snr_idx);
        y_ev1 = add_sensor_noise(y_ev1_clean, t, SNR_dB);
        y_ev2 = add_sensor_noise(y_ev2_clean, t, SNR_dB);

        T_noisy = table(t, y_ev1(:, 1), y_ev1(:, 2), y_ev2(:, 1), y_ev2(:, 2), ...
            'VariableNames', {'time', 'id_ev1', 'iq_ev1', 'id_ev2', 'iq_ev2'});
        fname = sprintf('sys_%05d_SNR_%d.csv', sys_id, SNR_dB);
        writetable(T_noisy, fullfile(out_dir, fname));
    end

    % Export ground-truth admittance and parameters
    ymat_fname = sprintf('sys_%05d_Ytrue.mat', sys_id);
    save(fullfile(out_dir, ymat_fname), 'Y_true', 'freqs_Hz', ...
         'Lf', 'Rf', 'Kp_i', 'Ki_i', 'Kp_pll', 'Ki_pll', 'Id_eq', 'Iq_eq');

    manifest(end+1, :) = {sys_id, Id_eq, Iq_eq, ymat_fname}; %#ok<SAGROW>

    if mod(sys_id, 200) == 0
        fprintf("  • Generated %4d / %4d systems (Stability success: %.1f%%)\n", ...
            sys_id, n_systems, (sys_id / attempts) * 100);
    end
end

if isempty(manifest)
    error("[Execution Failed] Manifest is empty. Verify state-space matrix formulation.");
end

%% 3. EXPORT MANIFEST FILE
manifest_table = cell2table(manifest, 'VariableNames', {'sys_id', 'Id_eq', 'Iq_eq', 'ytrue_file'});
writetable(manifest_table, fullfile(out_dir, "manifest.csv"));

fprintf("=========================================================================\n");
fprintf("Generation complete: %d systems written to '%s/'\n", sys_id, out_dir);
fprintf("SNR sweep levels: %s dB\n", mat2str(SNR_dB_list));
fprintf("=========================================================================\n");

%% ------------------------------------------------------------------------
% LOCAL FUNCTION: IN-SITU SENSOR NOISE INJECTION
% ------------------------------------------------------------------------
function y_noisy = add_sensor_noise(y_clean, t, SNR_dB)
    % Injects zero-mean Gaussian measurement noise scaled relative to 
    % the post-perturbation signal variance on each channel.
    y_noisy = y_clean;
    step_idx = (t >= 0.1);

    for ch = 1:size(y_clean, 2)
        sig = y_clean(step_idx, ch);
        sig_power = mean((sig - mean(sig)).^2);
        if sig_power == 0
            sig_power = eps;
        end
        noise_power = sig_power / (10^(SNR_dB / 10));
        noise = sqrt(noise_power) * randn(size(y_clean, 1), 1);
        y_noisy(:, ch) = y_clean(:, ch) + noise;
    end
end
