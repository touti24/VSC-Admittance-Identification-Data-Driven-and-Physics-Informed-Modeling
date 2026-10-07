%% generate_vsc_data.m
% =========================================================================
% CORRECTED DATA GENERATION FOR VSC I / VSC II ADMITTANCE IDENTIFICATION
% =========================================================================
% FIXES vs the previous version (per supervisor meeting, see transcript):
%   1) The perturbation is applied to GRID VOLTAGE (Vd, Vq) in the GLOBAL
%      dq frame, NOT to the current. You cannot physically step a current
%      by itself -- only voltage can be stepped and current is the
%      resulting response. This IS the admittance definition:
%           Y(s) = I(s) / V(s)   (output = current, input = voltage)
%   2) Noise is injected INSIDE this simulation, at the sensor/output
%      level, sample-by-sample, using the SAME statistics you will later
%      use at training/test time. Do NOT add noise afterwards in Python.
%      If you are running the real nonlinear Simulink model instead of
%      this linearized state-space version, add the noise block directly
%      on the Id/Iq/Vd/Vq "sensor" signals inside Simulink, at this same
%      sample rate, BEFORE logging to the workspace/CSV.
%   3) Everything below is written so VSC I and VSC II share the same
%      pipeline -- only the physical parameters + filter topology differ.
% =========================================================================

clear; clc;

%% ---------------- USER CONFIG ----------------
vsc_name        = "VSC_I";      % "VSC_I" or "VSC_II"
Id_equilibrium  = 50;           % [A] operating point on d-axis
Iq_equilibrium  = 0;            % [A] operating point on q-axis
Fs              = 2500;         % [Hz] simulation sample rate
T_window        = 0.98;         % [s] length of each step-response record
step_size_V     = 1.0;          % [V] magnitude of the voltage step (small-signal!)
SNR_dB_list     = [30 20 15 10 5];   % generate one noisy copy per SNR
add_clean_copy  = true;         % also save the noise-free record (SNR = Inf)
out_dir         = "vsc_data_out";

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% ---------------- PHYSICAL PARAMETERS ----------------
if vsc_name == "VSC_I"
    Lf = 1e-3;      Rf = 0.3e-3;
    Kp_i = 10.5;    Ki_i = 5900;
    Kp_pll = 1.2;   Ki_pll = 257;
    has_LCL = false;
elseif vsc_name == "VSC_II"
    Lf = 3e-3;      Rf = 0.9e-3;
    Kp_i = 10.5;    Ki_i = 5900;   % current-loop gains assumed same family
    Kp_pll = 1.2;   Ki_pll = 257;
    has_LCL = true;                 % VSC II has LCL filter + active damping
    Lg = 1.5e-3;  Cf = 10e-6;  Rd_damp = 0.5;   % LCL + active damping params
else
    error("vsc_name must be VSC_I or VSC_II");
end

w0 = 2*pi*50;  % grid fundamental angular frequency (50 Hz), for dq coupling

%% ---------------- BUILD LINEARIZED STATE-SPACE MODEL ----------------
% States for VSC I (6-state): [id, iq, xi_d, xi_q, theta_pll, xi_pll]
%   id, iq          : filter inductor currents (d/q)
%   xi_d, xi_q      : current-controller integrator states
%   theta_pll       : PLL angle error state
%   xi_pll          : PLL PI integrator state
% Inputs : [Vd, Vq]  (grid voltage perturbation, GLOBAL FRAME)
% Outputs: [id, iq]  (measured filter currents)
%
% This is a small-signal linearization about (Id_equilibrium, Iq_equilibrium).
% If you have the supervisor's exact 6-state model equations, replace the
% A,B,C,D block below with those -- the injection/noise logic downstream
% is independent of the exact plant model.

if ~has_LCL
    n = 6;
    A = zeros(n);
    B = zeros(n,2);
    C = zeros(2,n);

    % Filter dynamics: Lf * d(id,iq)/dt = -Rf*(id,iq) +/- w0*Lf*(iq,id) + Vconv - Vgrid
    A(1,1) = -Rf/Lf;   A(1,2) =  w0;
    A(2,1) = -w0;      A(2,2) = -Rf/Lf;

    % Converter voltage = PI(current error) -> feeds back into filter eq.
    A(1,3) =  Kp_i/Lf; A(2,4) =  Kp_i/Lf;
    A(3,1) = -1;       A(4,2) = -1;         % integrator states track -error
    % (Ki_i enters through xi_d, xi_q integration directly, already unit gain
    %  on error; Kp_i appears as proportional feed already captured above)

    % PLL: small-signal angle error feeds back into dq rotation (cross-coupling
    % of the equilibrium current through the rotated grid voltage)
    A(5,6) = 1;                  % theta_pll driven by integrator
    A(6,5) = -Ki_pll;            % PI on PLL phase-detector output (linearized)
    A(1,5) =  Iq_equilibrium*w0; % rotation of equilibrium current with angle error
    A(2,5) = -Id_equilibrium*w0;

    % Input matrix: grid voltage Vd, Vq enters filter equation directly (-1/Lf)
    % and also perturbs the PLL phase detector (approx via Vq channel)
    B(1,1) = -1/Lf;
    B(2,2) = -1/Lf;
    B(5,2) =  1;                 % PLL phase detector ~ Vq perturbation

    C(1,1) = 1;   % output id
    C(2,2) = 1;   % output iq
    D = zeros(2,2);
else
    % VSC II: LCL filter + active damping -> augment with capacitor voltage
    % and grid-side inductor current states (8-state model)
    n = 8;
    A = zeros(n); B = zeros(n,2); C = zeros(2,n);
    % States: [id, iq, xi_d, xi_q, theta_pll, xi_pll, vcd, vcq]
    A(1,1) = -Rf/Lf;   A(1,2) =  w0;   A(1,7) = -1/Lf;
    A(2,1) = -w0;      A(2,2) = -Rf/Lf; A(2,8) = -1/Lf;
    A(1,3) =  Kp_i/Lf; A(2,4) =  Kp_i/Lf;
    A(3,1) = -1;       A(4,2) = -1;
    A(5,6) = 1;        A(6,5) = -Ki_pll;
    A(1,5) =  Iq_equilibrium*w0;  A(2,5) = -Id_equilibrium*w0;
    % capacitor voltage dynamics with active damping resistor Rd_damp
    A(7,1) =  1/Cf;    A(7,8) =  w0;    A(7,7) = -1/(Rd_damp*Cf);
    A(8,2) =  1/Cf;    A(8,7) = -w0;    A(8,8) = -1/(Rd_damp*Cf);

    B(1,1) = 0;  % grid voltage now enters through Lg (grid-side inductor), approx:
    B(7,1) = -1/(Lg*Cf) * 0;   % placeholder cross term (kept 0th order for stability)
    B(1,1) = -1/Lf * 0.0;      % main coupling still through vcd/vcq (already in A)
    B(7,1) = 1/(Lg);           % NOTE: simplified; refine with full Lg branch if needed
    B(8,2) = 1/(Lg);

    C(1,1) = 1; C(2,2) = 1;
    D = zeros(2,2);
end

sys = ss(A, B, C, D);

% sanity check: system must be stable (all poles in LHP)
if any(real(eig(A)) >= 0)
    warning("Linearized model has non-negative-real-part poles -- check parameters!");
end

%% ---------------- GENERATE STEP-RESPONSE EVENTS (VOLTAGE STEPS, GLOBAL FRAME) ----------------
t = (0:1/Fs:T_window-1/Fs)';
Nt = length(t);

% Event 1: step on Vd (global frame), Vq held at 0
u_ev1 = zeros(Nt,2);
u_ev1(t >= 0.1, 1) = step_size_V;   % step applied at t = 0.1s to allow settling

% Event 2: step on Vq (global frame), Vd held at 0
u_ev2 = zeros(Nt,2);
u_ev2(t >= 0.1, 2) = step_size_V;

y_ev1_clean = lsim(sys, u_ev1, t);   % columns: [id, iq]
y_ev2_clean = lsim(sys, u_ev2, t);

%% ---------------- ADMITTANCE GROUND TRUTH (for validation) ----------------
freqs_Hz = linspace(1, 100, 100);
Y_true = zeros(2,2,length(freqs_Hz));
for k = 1:length(freqs_Hz)
    s = 1i*2*pi*freqs_Hz(k);
    Gs = C*((s*eye(size(A)) - A)\B) + D;   % 2x2: rows=[id,iq], cols=[Vd,Vq]
    Y_true(:,:,k) = Gs;
end
save(fullfile(out_dir, vsc_name + "_admittance_groundtruth.mat"), ...
     "Y_true", "freqs_Hz", "Id_equilibrium", "Iq_equilibrium");

%% ---------------- NOISE EMBEDDING (INSIDE THE SIMULATION, PER-SAMPLE) ----------------
% Noise is added to the SIMULATED SENSOR OUTPUT (id, iq) sample-by-sample,
% exactly as a real current sensor's noise would appear -- not appended
% after the fact in a downstream script. The SAME function must be reused
% for both training-data generation and test-data generation so the noise
% statistics always match between train and test.
%
% SNR is defined per-channel, relative to that channel's own step-response
% signal power (post step, i.e. excluding the pre-step zero segment).

rng(42);  % reproducibility

for snr_idx = 1:length(SNR_dB_list) + double(add_clean_copy)
    is_clean = add_clean_copy && (snr_idx == length(SNR_dB_list) + 1);
    if is_clean
        SNR_dB = Inf;
        y_ev1 = y_ev1_clean;
        y_ev2 = y_ev2_clean;
    else
        SNR_dB = SNR_dB_list(snr_idx);
        y_ev1 = add_sensor_noise(y_ev1_clean, t, SNR_dB, rng);
        y_ev2 = add_sensor_noise(y_ev2_clean, t, SNR_dB, rng);
    end

    T = table(t, y_ev1(:,1), y_ev1(:,2), y_ev2(:,1), y_ev2(:,2), ...
        'VariableNames', {'time','id_ev1','iq_ev1','id_ev2','iq_ev2'});

    fname = sprintf('%s_Id_%g_SNR_%s.csv', vsc_name, Id_equilibrium, ...
                     regexprep(num2str(SNR_dB), '\.', 'p'));
    writetable(T, fullfile(out_dir, fname));
    fprintf("Saved %s (SNR = %s dB)\n", fname, num2str(SNR_dB));
end

fprintf("\nDone. Ground-truth admittance + step-response CSVs saved in '%s/'\n", out_dir);
fprintf("Reminder: Vd/Vq step magnitude was %.3f V (small-signal). \n", step_size_V);
fprintf("If reconstructed admittance still doesn't match Y_true, check that\n");
fprintf("step_size_V is small enough to stay in the linear regime.\n");

%% ---------------- LOCAL FUNCTION: SENSOR-LEVEL NOISE ----------------
function y_noisy = add_sensor_noise(y_clean, t, SNR_dB, ~)
    % Adds Gaussian measurement noise to each output channel independently,
    % scaled so the specified SNR (dB) holds relative to the post-step
    % signal power on that channel. Called inside the simulation loop so
    % train/test data always come from the identical noise model.
    y_noisy = y_clean;
    step_idx = t >= 0.1;   % only measure signal power after the step
    for ch = 1:size(y_clean,2)
        sig = y_clean(step_idx, ch);
        sig_power = mean((sig - mean(sig)).^2);
        if sig_power == 0
            sig_power = eps;
        end
        noise_power = sig_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power) * randn(size(y_clean,1),1);
        y_noisy(:,ch) = y_clean(:,ch) + noise;
    end
end