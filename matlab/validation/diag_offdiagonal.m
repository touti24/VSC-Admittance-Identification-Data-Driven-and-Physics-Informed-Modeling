% =========================================================================
% diag_offdiagonal.m
%
%  Checks whether the cross-coupling voltage responses (Vq from q-injection,
%  Vd from d-injection) are genuine signal or sitting on the noise floor.
%  Run this BEFORE trusting Zdq / Zqd values.
%
%  Prints SNR at fp for each signal and overlays time-domain plots so you
%  can visually confirm whether a sine is present after the injection gate.
% =========================================================================

MODEL_NAME = 'VSC_Iavrg';
fp_test    = 20;    % [Hz]
Id_test    = 50;    % [A]
Iq_test    = 0;     % [A]

run('params_VSC_I.m');

T_SETTLE  = 8.0;
N_cycles  = 20;                        % more cycles = better frequency resolution
T_perturb = N_cycles / fp_test;
N_samples = round(T_perturb / Ts);
T_total   = T_SETTLE + T_perturb;

t_vector       = (0 : Ts : T_total)';
injection_gate = double(t_vector >= T_SETTLE);
amp            = 0.05 * Id_test;       % nominal 5%

load_system(MODEL_NAME);
assignin('base', 'Id_ref', Id_test);
assignin('base', 'Iq_ref', Iq_test);
simOpts = simset('SrcWorkspace', 'base');

fprintf('=========================================================\n');
fprintf('Off-diagonal SNR diagnostic at fp=%d Hz, Id=%d A\n\n', fp_test, Id_test);

% ── INJECTION 1: q-axis (paper d-axis) ───────────────────────────────────
id_sig = zeros(size(t_vector));
iq_sig = amp * sin(2*pi*fp_test * t_vector) .* injection_gate;

assignin('base', 'id_perturb', timeseries(id_sig, t_vector));
assignin('base', 'iq_perturb', timeseries(iq_sig, t_vector));

simOut1 = sim(MODEL_NAME, [0 T_total], simOpts);

t1    = simOut1.tout;
Iq1_t = extract_column(simOut1.Id_Iq_data, 2);   % injected channel
Vd1_t = extract_column(simOut1.Vd_Vq_data, 1);   % cross-coupling response
Vq1_t = extract_column(simOut1.Vd_Vq_data, 2);   % direct response (should be larger)

fprintf('--- Injection 1: q-axis sim (d-axis paper) ---\n');
snr_Vd1 = compute_snr(Vd1_t, t1, fp_test, T_SETTLE, N_samples, Ts);
snr_Vq1 = compute_snr(Vq1_t, t1, fp_test, T_SETTLE, N_samples, Ts);
fprintf('  Vd1 (cross)  SNR @ %d Hz = %.1f dB\n', fp_test, snr_Vd1);
fprintf('  Vq1 (direct) SNR @ %d Hz = %.1f dB\n\n', fp_test, snr_Vq1);

% ── INJECTION 2: d-axis (paper q-axis) ───────────────────────────────────
id_sig = amp * sin(2*pi*fp_test * t_vector) .* injection_gate;
iq_sig = zeros(size(t_vector));

assignin('base', 'id_perturb', timeseries(id_sig, t_vector));
assignin('base', 'iq_perturb', timeseries(iq_sig, t_vector));

simOut2 = sim(MODEL_NAME, [0 T_total], simOpts);

t2    = simOut2.tout;
Id2_t = extract_column(simOut2.Id_Iq_data, 1);   % injected channel
Vd2_t = extract_column(simOut2.Vd_Vq_data, 1);   % direct response
Vq2_t = extract_column(simOut2.Vd_Vq_data, 2);   % cross-coupling response

fprintf('--- Injection 2: d-axis sim (q-axis paper) ---\n');
snr_Vd2 = compute_snr(Vd2_t, t2, fp_test, T_SETTLE, N_samples, Ts);
snr_Vq2 = compute_snr(Vq2_t, t2, fp_test, T_SETTLE, N_samples, Ts);
fprintf('  Vd2 (direct) SNR @ %d Hz = %.1f dB\n', fp_test, snr_Vd2);
fprintf('  Vq2 (cross)  SNR @ %d Hz = %.1f dB\n\n', fp_test, snr_Vq2);

fprintf('=========================================================\n');
fprintf('Interpretation guide:\n');
fprintf('  SNR > 30 dB  → reliable phasor, trust the result\n');
fprintf('  SNR 15-30 dB → marginal, increase N_cycles or amplitude\n');
fprintf('  SNR < 15 dB  → noise dominated, off-diagonal unreliable\n\n');

% ── Time-domain plots ─────────────────────────────────────────────────────
i0_1 = find(t1 >= T_SETTLE, 1);
i0_2 = find(t2 >= T_SETTLE, 1);

figure(20); clf;
sgtitle(sprintf('Injection 1 (q-sim/d-paper) — fp=%d Hz  amp=%.2f A', fp_test, amp));

subplot(2,2,1);
plot(t1, Iq1_t); xline(T_SETTLE,'r--');
title('Iq1 — injected (should show clear sine)'); grid on; ylabel('A');

subplot(2,2,2);
plot(t1, Vq1_t); xline(T_SETTLE,'r--');
title('Vq1 — direct voltage response'); grid on; ylabel('V');

subplot(2,2,3);
plot(t1(i0_1:end), Vd1_t(i0_1:end)); 
title(sprintf('Vd1 — cross-coupling (SNR=%.1f dB)', snr_Vd1)); grid on; ylabel('V');
xlabel('time (s)');

subplot(2,2,4);
% FFT of Vd1 steady-state segment
seg = Vd1_t(i0_1 : i0_1+N_samples-1);
seg = (seg - mean(seg)) .* hann(length(seg));
NFFT = 2^nextpow2(length(seg)*4);
F    = abs(fft(seg, NFFT)) * 2/sum(hann(N_samples));
fs   = 1/Ts;
f_ax = (0:NFFT-1)*fs/NFFT;
plot(f_ax(1:NFFT/2), F(1:NFFT/2));
xline(fp_test,'r--',sprintf('%d Hz',fp_test));
xlim([0 200]); grid on;
title('FFT of Vd1 (cross) — look for spike at fp'); ylabel('V'); xlabel('Hz');

figure(21); clf;
sgtitle(sprintf('Injection 2 (d-sim/q-paper) — fp=%d Hz  amp=%.2f A', fp_test, amp));

subplot(2,2,1);
plot(t2, Id2_t); xline(T_SETTLE,'r--');
title('Id2 — injected (should show clear sine)'); grid on; ylabel('A');

subplot(2,2,2);
plot(t2, Vd2_t); xline(T_SETTLE,'r--');
title('Vd2 — direct voltage response'); grid on; ylabel('V');

subplot(2,2,3);
plot(t2(i0_2:end), Vq2_t(i0_2:end));
title(sprintf('Vq2 — cross-coupling (SNR=%.1f dB)', snr_Vq2)); grid on; ylabel('V');
xlabel('time (s)');

subplot(2,2,4);
seg = Vq2_t(i0_2 : i0_2+N_samples-1);
seg = (seg - mean(seg)) .* hann(length(seg));
NFFT = 2^nextpow2(length(seg)*4);
F    = abs(fft(seg, NFFT)) * 2/sum(hann(N_samples));
f_ax = (0:NFFT-1)*fs/NFFT;
plot(f_ax(1:NFFT/2), F(1:NFFT/2));
xline(fp_test,'r--',sprintf('%d Hz',fp_test));
xlim([0 200]); grid on;
title('FFT of Vq2 (cross) — look for spike at fp'); ylabel('V'); xlabel('Hz');

close_system(MODEL_NAME, 0);


% =========================================================================
% Compute SNR at fp relative to surrounding noise floor
%   SNR = 20*log10( amplitude_at_fp / rms_noise_floor )
% =========================================================================
function snr_dB = compute_snr(sig, t, fp, T_settle, N_samples, Ts)
    i0  = find(t >= T_settle, 1);
    seg = sig(i0 : i0 + N_samples - 1);
    seg = (seg - mean(seg)) .* hann(N_samples);

    NFFT  = 2^nextpow2(N_samples * 4);
    fs    = 1 / Ts;
    f_ax  = (0 : NFFT-1) * fs / NFFT;
    F_mag = abs(fft(seg, NFFT));

    [~, k_fp] = min(abs(f_ax - fp));
    sig_amp   = F_mag(k_fp);

    % Noise floor: median of bins excluding ±5 bins around fp
    excl  = max(1, k_fp-5) : min(NFFT, k_fp+5);
    noise_bins        = F_mag(1 : NFFT/2);
    noise_bins(excl)  = NaN;
    noise_floor       = median(noise_bins, 'omitnan');

    if noise_floor < eps
        snr_dB = Inf;
    else
        snr_dB = 20 * log10(sig_amp / noise_floor);
    end
end


% =========================================================================
% HELPER: extract one column from To-Workspace output
% =========================================================================
function col = extract_column(sig, idx)
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