%% generate_vsc_data.m
% =========================================================================
% Small-Signal Admittance Dataset Generator for Grid-Connected VSCs
%
% Generates time-domain step responses and analytical ground-truth admittance
% for VSC I (and VSC II). Preserves exact numerical compatibility with
% Notebook 02 and Notebook 03 benchmarks.
% =========================================================================

clear; clc; close all;

%% 1. USER CONFIGURATION
vsc_name        = "VSC_I";            % "VSC_I" or "VSC_II"
Id_equilibrium  = 50;                 % [A] Operating point on d-axis
Iq_equilibrium  = 0;                  % [A] Operating point on q-axis
Fs              = 2500;               % [Hz] Simulation sample rate
T_window        = 0.98;               % [s] Length of each step-response record
step_size_V     = 1.0;                % [V] Small-signal voltage step perturbation
SNR_dB_list     = [30, 20, 15, 10, 5];% Target SNR levels [dB]
add_clean_copy  = true;               % Save noise-free record (SNR = Inf)
out_dir         = "vsc_data_out";     % Output directory

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% 2. PHYSICAL PARAMETERS
w0    = 2 * pi * 50;                  % Fundamental grid frequency (50 Hz)
Vd_ss = 311.0;                        % Nominal peak phase voltage (V)

if vsc_name == "VSC_I"
    Lf      = 1.0e-3;   Rf      = 0.3e-3;   % Filter inductance [H] & resistance [Ohm]
    Kp_i    = 10.5;     Ki_i    = 5900.0;   % Current loop PI gains
    Kp_pll  = 1.2;      Ki_pll  = 257.0;    % PLL PI gains
    has_LCL = false;
elseif vsc_name == "VSC_II"
    Lf      = 3.0e-3;   Rf      = 0.9e-3;
    Kp_i    = 10.5;     Ki_i    = 5900.0;
    Kp_pll  = 1.2;      Ki_pll  = 257.0;
    has_LCL = true;
    Lg      = 1.5e-3;   Cf      = 10.0e-6;  Rd_damp = 0.5;
else
    error("vsc_name must be 'VSC_I' or 'VSC_II'");
end

%% 3. LINEARIZED STATE-SPACE MODEL
if ~has_LCL
    % VSC I (6-state model matching training_data_out5 and Notebook 02)
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
    A(1, 5) =  Iq_equilibrium * w0;
    A(2, 5) = -Id_equilibrium * w0;

    B(1, 1) = -1.0 / Lf;
    B(2, 2) = -1.0 / Lf;
    B(5, 2) =  1.0;

    C(1, 1) = 1.0;
    C(2, 2) = 1.0;
else
    % VSC II (8-state LCL model with active damping)
    n = 8;
    A = zeros(n, n);
    B = zeros(n, 2);
    C = zeros(2, n);
    D = zeros(2, 2);

    A(1, 1) = -Rf / Lf;      A(1, 2) =  w0;          A(1, 7) = -1.0 / Lf;
    A(2, 1) = -w0;           A(2, 2) = -Rf / Lf;      A(2, 8) = -1.0 / Lf;
    A(1, 3) =  Kp_i / Lf;    A(2, 4) =  Kp_i / Lf;
    A(3, 1) = -1.0;          A(4, 2) = -1.0;

    A(5, 5) = -Kp_pll * Vd_ss;
    A(5, 6) =  1.0;
    A(6, 5) = -Ki_pll * Vd_ss;
    A(1, 5) =  Iq_equilibrium * w0;
    A(2, 5) = -Id_equilibrium * w0;

    A(7, 1) =  1.0 / Cf;     A(7, 8) =  w0;          A(7, 7) = -1.0 / (Rd_damp * Cf);
    A(8, 2) =  1.0 / Cf;     A(8, 7) = -w0;          A(8, 8) = -1.0 / (Rd_damp * Cf);

    B(5, 2) =  1.0;
    B(7, 1) =  1.0 / Lg;
    B(8, 2) =  1.0 / Lg;

    C(1, 1) = 1.0;
    C(2, 2) = 1.0;
end

sys = ss(A, B, C, D);

if any(real(eig(A)) >= 0)
    warning("Linearized model has non-negative real part poles.");
end

%% 4. SIMULATE VOLTAGE-STEP TRANSIENTS (TWO-EVENT PROTOCOL)
dt = 1.0 / Fs;
t  = (0 : dt : T_window - dt)';
Nt = length(t);

% Event 1: Step perturbation on Vd at t >= 0.1 s
u_ev1 = zeros(Nt, 2);
u_ev1(t >= 0.1, 1) = step_size_V;

% Event 2: Step perturbation on Vq at t >= 0.1 s
u_ev2 = zeros(Nt, 2);
u_ev2(t >= 0.1, 2) = step_size_V;

y_ev1_clean = lsim(sys, u_ev1, t);   % Outputs: [id_ev1, iq_ev1]
y_ev2_clean = lsim(sys, u_ev2, t);   % Outputs: [id_ev2, iq_ev2]

%% 5. ANALYTICAL GROUND-TRUTH ADMITTANCE (EXACT 100-POINT GRID)
% PRESERVED EXACTLY: 100 linearly spaced points from 1 to 100 Hz
freqs_Hz = linspace(1, 100, 100);
Y_true   = zeros(2, 2, length(freqs_Hz));

for k = 1:length(freqs_Hz)
    s = 1j * 2 * pi * freqs_Hz(k);
    Y_true(:, :, k) = C * ((s * eye(n) - A) \ B) + D;
end

gt_file = fullfile(out_dir, vsc_name + "_admittance_groundtruth.mat");
save(gt_file, "Y_true", "freqs_Hz", "Id_equilibrium", "Iq_equilibrium", ...
     "Lf", "Rf", "Kp_i", "Ki_i", "Kp_pll", "Ki_pll");
fprintf("[Output] Saved: %s\n", gt_file);

%% 6. NOISE EMBEDDING & CSV EXPORT (IDENTICAL SEED & CALL SEQUENCE)
rng(42);  % Preserves identical pseudo-random sequence

for snr_idx = 1 : (length(SNR_dB_list) + double(add_clean_copy))
    is_clean = add_clean_copy && (snr_idx == length(SNR_dB_list) + 1);

    if is_clean
        SNR_dB  = Inf;
        y_ev1   = y_ev1_clean;
        y_ev2   = y_ev2_clean;
        snr_str = "Inf";
    else
        SNR_dB  = SNR_dB_list(snr_idx);
        y_ev1   = add_sensor_noise(y_ev1_clean, t, SNR_dB);
        y_ev2   = add_sensor_noise(y_ev2_clean, t, SNR_dB);
        snr_str = num2str(SNR_dB);
    end

    T = table(t, y_ev1(:, 1), y_ev1(:, 2), y_ev2(:, 1), y_ev2(:, 2), ...
        'VariableNames', {'time', 'id_ev1', 'iq_ev1', 'id_ev2', 'iq_ev2'});

    fname = sprintf('%s_Id_%g_SNR_%s.csv', vsc_name, Id_equilibrium, snr_str);
    writetable(T, fullfile(out_dir, fname));
    fprintf("Saved %s (SNR = %s dB)\n", fname, snr_str);
end

fprintf("\nDone. Datasets written to '%s/'.\n", out_dir);

%% ------------------------------------------------------------------------
% LOCAL FUNCTION: SENSOR-LEVEL NOISE (EXACT ORIGINAL LOGIC)
% ------------------------------------------------------------------------
function y_noisy = add_sensor_noise(y_clean, t, SNR_dB)
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
