%% sanity_check_amplitude_adapted.m
% =========================================================================
% Small-Signal Perturbation Amplitude Linearity Verification
%
% Description:
%   Validates small-signal linearity for the average-value VSC model
%   (VSC_Iavrg) across three perturbation injection amplitudes:
%     - Half     : 2.5% of nominal current (Id = 50 A)
%     - Nominal  : 5.0% of nominal current
%     - Double   : 10.0% of nominal current
%
% Methodology:
%   - Excitation: Injects sinusoidal perturbations at test frequency fp = 20 Hz
%     across both orthogonal channels (q-sim/d-paper and d-sim/q-paper).
%   - Extraction: Extracts complex phasors using Hann-windowed FFTs.
%   - Matrix Solving: Computes the 2x2 impedance matrix Z(s) = inv(Y(s)).
%   - Verification: Checks that impedance variations between half and double
%     amplitudes remain within 5% to guarantee operation in the linear regime.
% =========================================================================

clear; clc; close all;

%% 1. OPERATING POINT & TEST CONFIGURATION
MODEL_NAME = 'VSC_Iavrg';
fp_test    = 20;                              % Test perturbation frequency [Hz]
Id_test    = 50;                              % Active current setpoint [A]
Iq_test    = 0;                               % Reactive current setpoint [A]

run('params_VSC_I.m');

T_SETTLE   = 8.0;                             % Settling time before perturbation [s]
load_system(MODEL_NAME);
assignin('base', 'Id_ref', Id_test);
assignin('base', 'Iq_ref', Iq_test);

% Perturbation amplitudes (2.5%, 5.0%, and 10.0% of steady-state active current)
amp_nominal = 0.05 * Id_test;
amp_half    = amp_nominal / 2;
amp_double  = amp_nominal * 2;

amps   = [amp_half, amp_nominal, amp_double];
labels = {'Half (2.5%)', 'Nominal (5%)', 'Double (10%)'};

fprintf('=========================================================\n');
fprintf('Linearity Sanity Check at fp = %d Hz, Id = %d A\n', fp_test, Id_test);
fprintf('Amplitudes: %.3f A (2.5%%), %.3f A (5%%), %.3f A (10%%)\n', ...
        amp_half, amp_nominal, amp_double);
fprintf('=========================================================\n\n');

%% 2. SIMULATION TIMING & BUFFER ALLOCATION
simOpts   = simset('SrcWorkspace', 'base');
N_cycles  = 5;                                % Evaluation window cycles
T_perturb = N_cycles / fp_test;
N_samples = round(T_perturb / Ts);
T_total   = T_SETTLE + T_perturb;

t_vector       = (0 : Ts : T_total)';
injection_gate = double(t_vector >= T_SETTLE);

% Pre-allocate impedance results across test amplitudes (Rows: Zdd, Zdq, Zqd, Zqq)
Z_results = zeros(4, 3);

%% 3. MULTI-AMPLITUDE PERTURBATION SWEEP
for k = 1:3
    fprintf('---------------------------------------------------------\n');
    fprintf('Case %d: %s (Amplitude = %.3f A)\n', k, labels{k}, amps(k));

    % ── Injection 1: q-sim = d-paper ─────────────────────────────────────
    assignin('base', 'id_perturb', timeseries(zeros(size(t_vector)), t_vector));
    assignin('base', 'iq_perturb', timeseries(amps(k) * sin(2*pi*fp_test*t_vector) .* injection_gate, t_vector));

    simOut1 = sim(MODEL_NAME, [0, T_total], simOpts);
    [Vd1, Vq1, Id1, Iq1] = extract_phasors_local(simOut1, fp_test, Ts, T_SETTLE, N_samples);

    % Optional diagnostic inspection on initial test amplitude
    if k == 1
        plot_diagnostic(simOut1, T_SETTLE, 10, 'Injection 1 (q-sim / d-paper)');
        fprintf('  DC Check: Id = %.3f A | Iq = %.3f A | Vd = %.3f V | Vq = %.3f V\n', ...
            mean(extract_column(simOut1.Id_Iq_data, 1)), ...
            mean(extract_column(simOut1.Id_Iq_data, 2)), ...
            mean(extract_column(simOut1.Vd_Vq_data, 1)), ...
            mean(extract_column(simOut1.Vd_Vq_data, 2)));
    end

    % ── Injection 2: d-sim = q-paper ─────────────────────────────────────
    assignin('base', 'id_perturb', timeseries(amps(k) * sin(2*pi*fp_test*t_vector) .* injection_gate, t_vector));
    assignin('base', 'iq_perturb', timeseries(zeros(size(t_vector)), t_vector));

    simOut2 = sim(MODEL_NAME, [0, T_total], simOpts);
    [Vd2, Vq2, Id2, Iq2] = extract_phasors_local(simOut2, fp_test, Ts, T_SETTLE, N_samples);

    % Display extracted phasors
    fprintf('  Phasors:\n');
    fprintf('    Iq1 = %+.4f %+.4fj | Vq1 = %+.4f %+.4fj | Vd1 = %+.4f %+.4fj\n', ...
        real(Iq1), imag(Iq1), real(Vq1), imag(Vq1), real(Vd1), imag(Vd1));
    fprintf('    Id2 = %+.4f %+.4fj | Vd2 = %+.4f %+.4fj | Vq2 = %+.4f %+.4fj\n', ...
        real(Id2), imag(Id2), real(Vd2), imag(Vd2), real(Vq2), imag(Vq2));

    % ── Matrix Assembly & Impedance Calculation ──────────────────────────
    % Target reference frame: Row 1 = d-paper (q-sim), Row 2 = q-paper (d-sim)
    I_mat = [Iq1, Iq2;
             Id1, Id2];

    V_mat = [Vq1, Vq2;
             Vd1, Vd2];

    rc = rcond(V_mat);
    fprintf('  Matrix Condition Number rcond(V_mat) = %.3e\n', rc);
    if rc < 1e-6
        warning('Voltage matrix is ill-conditioned (rcond < 1e-6). Skipping Case %d.', k);
        Z_results(:, k) = NaN;
        clear simOut1 simOut2;
        continue;
    end

    Y_sim = V_mat \ I_mat;                   % Admittance matrix
    Z_sim = inv(Y_sim);                      % Impedance matrix

    Zdd = Z_sim(1, 1);
    Zdq = Z_sim(1, 2);
    Zqd = Z_sim(2, 1);
    Zqq = Z_sim(2, 2);

    fprintf('  Impedance Response (dB):\n');
    fprintf('    Zdd: %.4f Ohm | Phase: %+.1f deg | %+.2f dB\n', abs(Zdd), angle(Zdd)*180/pi, 20*log10(abs(Zdd)));
    fprintf('    Zdq: %.4f Ohm | Phase: %+.1f deg | %+.2f dB\n', abs(Zdq), angle(Zdq)*180/pi, 20*log10(abs(Zdq)));
    fprintf('    Zqd: %.4f Ohm | Phase: %+.1f deg | %+.2f dB\n', abs(Zqd), angle(Zqd)*180/pi, 20*log10(abs(Zqd)));
    fprintf('    Zqq: %.4f Ohm | Phase: %+.1f deg | %+.2f dB\n', abs(Zqq), angle(Zqq)*180/pi, 20*log10(abs(Zqq)));

    Z_results(:, k) = [Zdd; Zdq; Zqd; Zqq];
    clear simOut1 simOut2;
end

%% 4. LINEARITY CONSISTENCY EVALUATION (PASS / FAIL)
fprintf('\n=========================================================\n');
fprintf('Linearity Consistency Check (Half vs. Double Amplitude):\n');

diff_hd = abs(abs(Z_results(:, 1)) - abs(Z_results(:, 3))) ./ ...
          (abs(Z_results(:, 2)) + eps) * 100;

names_Z = {'Zdd', 'Zdq', 'Zqd', 'Zqq'};
for i = 1:4
    fprintf('  %s Variation: %.2f%%\n', names_Z{i}, diff_hd(i));
end

if max(diff_hd) < 5
    fprintf('\nRESULT: PASS (Maximum variation = %.2f%% < 5%% threshold).\n', max(diff_hd));
    fprintf('Small-signal linearity assumption is strictly satisfied.\n');
else
    fprintf('\nRESULT: FAIL (Maximum variation = %.2f%% >= 5%% threshold).\n', max(diff_hd));
    fprintf('Nonlinear distortion detected. Consider lowering perturbation amplitude.\n');
end
fprintf('=========================================================\n');

close_system(MODEL_NAME, 0);

%% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function plot_diagnostic(simOut, T_SETTLE, fig_num, ttl)
    % Renders time-domain traces of current and voltage states
    Id_t = extract_column(simOut.Id_Iq_data, 1);
    Iq_t = extract_column(simOut.Id_Iq_data, 2);
    Vd_t = extract_column(simOut.Vd_Vq_data, 1);
    Vq_t = extract_column(simOut.Vd_Vq_data, 2);
    t    = simOut.tout;

    figure(fig_num); clf;
    subplot(2, 2, 1); plot(t, Id_t); xline(T_SETTLE, 'r--'); title('Id'); grid on; ylabel('A');
    subplot(2, 2, 2); plot(t, Iq_t); xline(T_SETTLE, 'r--'); title('Iq'); grid on; ylabel('A');
    subplot(2, 2, 3); plot(t, Vd_t); xline(T_SETTLE, 'r--'); title('Vd'); grid on; ylabel('V'); xlabel('Time (s)');
    subplot(2, 2, 4); plot(t, Vq_t); xline(T_SETTLE, 'r--'); title('Vq'); grid on; ylabel('V'); xlabel('Time (s)');
    sgtitle(ttl);
end

function col = extract_column(sig, idx)
    % Extracts single column vector from timeseries or array structure
    if isa(sig, 'tscollection')
        names = sig.TimeSeriesNames;
        col   = sig.(names{idx}).Data;
    elseif isa(sig, 'timeseries')
        if size(sig.Data, 2) >= idx
            col = sig.Data(:, idx);
        else
            col = sig.Data(:, 1);
        end
    else
        col = sig(:, idx);
    end
    col = double(col(:));
end

function [Vd_ph, Vq_ph, Id_ph, Iq_ph] = extract_phasors_local(simOut, fp, Ts, T_settle, N_samples)
    % Extracts fundamental complex phasors at frequency fp via Hann-windowed FFT
    t      = simOut.tout;
    Id_all = extract_column(simOut.Id_Iq_data, 1);
    Iq_all = extract_column(simOut.Id_Iq_data, 2);
    Vd_all = extract_column(simOut.Vd_Vq_data, 1);
    Vq_all = extract_column(simOut.Vd_Vq_data, 2);

    i_start = find(t >= T_settle, 1, 'first');
    i_end   = i_start + N_samples - 1;
    if i_end > length(t)
        warning('Truncating extraction window to match simulation record.');
        i_end     = length(t);
        N_samples = i_end - i_start + 1;
    end

    w     = hann(N_samples);
    scale = 2 / sum(w);
    NFFT  = 2^nextpow2(N_samples * 4);
    fs    = 1 / Ts;
    f_ax  = (0 : NFFT - 1) * fs / NFFT;
    [~, k_fp] = min(abs(f_ax - fp));

    Vd_fft = fft((Vd_all(i_start:i_end) - mean(Vd_all(i_start:i_end))) .* w, NFFT);
    Vq_fft = fft((Vq_all(i_start:i_end) - mean(Vq_all(i_start:i_end))) .* w, NFFT);
    Id_fft = fft((Id_all(i_start:i_end) - mean(Id_all(i_start:i_end))) .* w, NFFT);
    Iq_fft = fft((Iq_all(i_start:i_end) - mean(Iq_all(i_start:i_end))) .* w, NFFT);

    Vd_ph = Vd_fft(k_fp) * scale;
    Vq_ph = Vq_fft(k_fp) * scale;
    Id_ph = Id_fft(k_fp) * scale;
    Iq_ph = Iq_fft(k_fp) * scale;
end
