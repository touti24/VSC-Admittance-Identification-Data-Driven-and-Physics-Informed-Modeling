%% diag_offdiagonal.m
% =========================================================================
% Cross-Coupling Signal-to-Noise Ratio (SNR) Diagnostic Script
%
% Description:
%   Evaluates whether the cross-coupling voltage responses (Vd from q-injection,
%   Vq from d-injection) represent genuine linear small-signal dynamics or
%   remain buried in the simulation/numerical noise floor.
%
% Verification Protocol:
%   1. Injects a 5% sinusoidal current perturbation at test frequency fp = 20 Hz
%      after steady-state settling (T_settle = 8.0 s).
%   2. Simulates orthogonal injections:
%        - Injection 1: q-sim (d-paper) perturbation -> measures Vd1 (cross) & Vq1 (direct)
%        - Injection 2: d-sim (q-paper) perturbation -> measures Vd2 (direct) & Vq2 (cross)
%   3. Computes spectral SNR at fp relative to the median noise floor:
%        SNR = 20 * log10( amplitude_at_fp / rms_noise_floor )
%   4. Renders time-domain traces and steady-state FFT spectra for visual confirmation.
% =========================================================================

clear; clc; close all;

%% 1. OPERATING POINT & TEST CONFIGURATION
MODEL_NAME = 'VSC_Iavrg';
fp_test    = 20;                              % Diagnostic perturbation frequency [Hz]
Id_test    = 50;                              % Active current setpoint [A]
Iq_test    = 0;                               % Reactive current setpoint [A]

run('params_VSC_I.m');

T_SETTLE   = 8.0;                             % System settling time [s]
N_cycles   = 20;                              % Perturbation cycles for spectral resolution
T_perturb  = N_cycles / fp_test;
N_samples  = round(T_perturb / Ts);
T_total    = T_SETTLE + T_perturb;

t_vector       = (0 : Ts : T_total)';
injection_gate = double(t_vector >= T_SETTLE);
amp            = 0.05 * Id_test;              % 5% nominal current amplitude [A]

load_system(MODEL_NAME);
assignin('base', 'Id_ref', Id_test);
assignin('base', 'Iq_ref', Iq_test);
simOpts = simset('SrcWorkspace', 'base');

fprintf('=========================================================\n');
fprintf('Off-Diagonal SNR Diagnostic: fp = %d Hz, Id = %d A\n', fp_test, Id_test);
fprintf('=========================================================\n\n');

%% 2. INJECTION 1: q-axis sim (d-axis paper)
id_sig = zeros(size(t_vector));
iq_sig = amp * sin(2 * pi * fp_test * t_vector) .* injection_gate;

assignin('base', 'id_perturb', timeseries(id_sig, t_vector));
assignin('base', 'iq_perturb', timeseries(iq_sig, t_vector));

simOut1 = sim(MODEL_NAME, [0, T_total], simOpts);

t1    = simOut1.tout;
Iq1_t = extract_column(simOut1.Id_Iq_data, 2);   % Injected current
Vd1_t = extract_column(simOut1.Vd_Vq_data, 1);   % Cross-coupling voltage response
Vq1_t = extract_column(simOut1.Vd_Vq_data, 2);   % Direct voltage response

snr_Vd1 = compute_snr(Vd1_t, t1, fp_test, T_SETTLE, N_samples, Ts);
snr_Vq1 = compute_snr(Vq1_t, t1, fp_test, T_SETTLE, N_samples, Ts);

fprintf('--- Injection 1: q-sim (d-paper) ---\n');
fprintf('  Vd1 (Cross-Coupling) SNR @ %d Hz = %.1f dB\n', fp_test, snr_Vd1);
fprintf('  Vq1 (Direct Channel) SNR @ %d Hz = %.1f dB\n\n', fp_test, snr_Vq1);

%% 3. INJECTION 2: d-axis sim (q-axis paper)
id_sig = amp * sin(2 * pi * fp_test * t_vector) .* injection_gate;
iq_sig = zeros(size(t_vector));

assignin('base', 'id_perturb', timeseries(id_sig, t_vector));
assignin('base', 'iq_perturb', timeseries(iq_sig, t_vector));

simOut2 = sim(MODEL_NAME, [0, T_total], simOpts);

t2    = simOut2.tout;
Id2_t = extract_column(simOut2.Id_Iq_data, 1);   % Injected current
Vd2_t = extract_column(simOut2.Vd_Vq_data, 1);   % Direct voltage response
Vq2_t = extract_column(simOut2.Vd_Vq_data, 2);   % Cross-coupling voltage response

snr_Vd2 = compute_snr(Vd2_t, t2, fp_test, T_SETTLE, N_samples, Ts);
snr_Vq2 = compute_snr(Vq2_t, t2, fp_test, T_SETTLE, N_samples, Ts);

fprintf('--- Injection 2: d-sim (q-paper) ---\n');
fprintf('  Vd2 (Direct Channel) SNR @ %d Hz = %.1f dB\n', fp_test, snr_Vd2);
fprintf('  Vq2 (Cross-Coupling) SNR @ %d Hz = %.1f dB\n\n', fp_test, snr_Vq2);

%% 4. DIAGNOSTIC INTERPRETATION SUMMARY
fprintf('=========================================================\n');
fprintf('SNR Interpretation Guidelines:\n');
fprintf('  • SNR > 30 dB  : High-precision phasor extraction (reliable).\n');
fprintf('  • SNR 15-30 dB : Acceptable; consider extending observation window.\n');
fprintf('  • SNR < 15 dB  : Noise-dominated; off-diagonal coupling unreliable.\n');
fprintf('=========================================================\n\n');

%% 5. TIME-DOMAIN & SPECTRAL VISUALIZATION
i0_1 = find(t1 >= T_SETTLE, 1);
i0_2 = find(t2 >= T_SETTLE, 1);

% ── Figure 20: Injection 1 Analysis ─────────────────────────────────────
figure(20); clf;
sgtitle(sprintf('Injection 1 (q-sim / d-paper) | fp = %d Hz, Amplitude = %.2f A', fp_test, amp));

subplot(2, 2, 1);
plot(t1, Iq1_t); xline(T_SETTLE, 'r--');
title('Iq1 - Injected Current'); grid on; ylabel('A');

subplot(2, 2, 2);
plot(t1, Vq1_t); xline(T_SETTLE, 'r--');
title('Vq1 - Direct Voltage Response'); grid on; ylabel('V');

subplot(2, 2, 3);
plot(t1(i0_1:end), Vd1_t(i0_1:end));
title(sprintf('Vd1 - Cross-Coupling (SNR = %.1f dB)', snr_Vd1));
grid on; ylabel('V'); xlabel('Time (s)');

subplot(2, 2, 4);
seg1 = Vd1_t(i0_1 : i0_1 + N_samples - 1);
seg1 = (seg1 - mean(seg1)) .* hann(length(seg1));
NFFT1 = 2^nextpow2(length(seg1) * 4);
F1    = abs(fft(seg1, NFFT1)) * 2 / sum(hann(N_samples));
fs    = 1 / Ts;
f_ax1 = (0 : NFFT1 - 1) * fs / NFFT1;
plot(f_ax1(1 : NFFT1 / 2), F1(1 : NFFT1 / 2));
xline(fp_test, 'r--', sprintf('%d Hz', fp_test));
xlim([0, 200]); grid on;
title('FFT Spectrum: Vd1 (Cross Channel)'); ylabel('V'); xlabel('Frequency (Hz)');

% ── Figure 21: Injection 2 Analysis ─────────────────────────────────────
figure(21); clf;
sgtitle(sprintf('Injection 2 (d-sim / q-paper) | fp = %d Hz, Amplitude = %.2f A', fp_test, amp));

subplot(2, 2, 1);
plot(t2, Id2_t); xline(T_SETTLE, 'r--');
title('Id2 - Injected Current'); grid on; ylabel('A');

subplot(2, 2, 2);
plot(t2, Vd2_t); xline(T_SETTLE, 'r--');
title('Vd2 - Direct Voltage Response'); grid on; ylabel('V');

subplot(2, 2, 3);
plot(t2(i0_2:end), Vq2_t(i0_2:end));
title(sprintf('Vq2 - Cross-Coupling (SNR = %.1f dB)', snr_Vq2));
grid on; ylabel('V'); xlabel('Time (s)');

subplot(2, 2, 4);
seg2 = Vq2_t(i0_2 : i0_2 + N_samples - 1);
seg2 = (seg2 - mean(seg2)) .* hann(length(seg2));
NFFT2 = 2^nextpow2(length(seg2) * 4);
F2    = abs(fft(seg2, NFFT2)) * 2 / sum(hann(N_samples));
f_ax2 = (0 : NFFT2 - 1) * fs / NFFT2;
plot(f_ax2(1 : NFFT2 / 2), F2(1 : NFFT2 / 2));
xline(fp_test, 'r--', sprintf('%d Hz', fp_test));
xlim([0, 200]); grid on;
title('FFT Spectrum: Vq2 (Cross Channel)'); ylabel('V'); xlabel('Frequency (Hz)');

close_system(MODEL_NAME, 0);

%% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function snr_dB = compute_snr(sig, t, fp, T_settle, N_samples, Ts)
    % Evaluates spectral signal-to-noise ratio at fp using a Hann window
    i0  = find(t >= T_settle, 1);
    seg = sig(i0 : i0 + N_samples - 1);
    seg = (seg - mean(seg)) .* hann(N_samples);

    NFFT  = 2^nextpow2(N_samples * 4);
    fs    = 1 / Ts;
    f_ax  = (0 : NFFT - 1) * fs / NFFT;
    F_mag = abs(fft(seg, NFFT));

    [~, k_fp] = min(abs(f_ax - fp));
    sig_amp   = F_mag(k_fp);

    % Exclude +/- 5 frequency bins around fp to evaluate noise floor
    excl = max(1, k_fp - 5) : min(NFFT, k_fp + 5);
    noise_bins = F_mag(1 : NFFT / 2);
    noise_bins(excl) = NaN;
    noise_floor = median(noise_bins, 'omitnan');

    if noise_floor < eps
        snr_dB = Inf;
    else
        snr_dB = 20 * log10(sig_amp / noise_floor);
    end
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
