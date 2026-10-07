% =========================================================================
% sanity_check_amplitude_adapted.m
%
%  Convention aligned to paper Fig.9:
%   - Paper d-axis = simulation q-axis
%   - Paper q-axis = simulation d-axis
%   - Both current and voltage rows swapped accordingly
%   - Plots IMPEDANCE Z in dB  (Fig.9 is Z not Y)
%   - Matrix formula: Y = V_mat \ I_mat,  Z = inv(Y)
% =========================================================================

MODEL_NAME = 'VSC_Iavrg';
fp_test    = 20;    % [Hz]
Id_test    = 50;    % [A]
Iq_test    = 0;     % [A]

run('params_VSC_I.m');

T_SETTLE  = 8.0;
load_system(MODEL_NAME);
assignin('base', 'Id_ref', Id_test);
assignin('base', 'Iq_ref', Iq_test);

amp_nominal = 0.05 * Id_test;
amp_half    = amp_nominal / 2;
amp_double  = amp_nominal * 2;

fprintf('=========================================================\n');
fprintf('Sanity check at fp=%d Hz, Id=%d A\n', fp_test, Id_test);
fprintf('Amplitudes: %.3f (half)  %.3f (nominal)  %.3f (double)\n\n', ...
        amp_half, amp_nominal, amp_double);

amps   = [amp_half, amp_nominal, amp_double];
labels = {'Half (2.5%)', 'Nominal (5%)', 'Double (10%)'};

% Rows: Zdd, Zdq, Zqd, Zqq  (paper frame)
Z_results = zeros(4, 3);

simOpts   = simset('SrcWorkspace', 'base');
N_cycles  = 5;
T_perturb = N_cycles / fp_test;
N_samples = round(T_perturb / Ts);
T_total   = T_SETTLE + T_perturb;

t_vector       = (0 : Ts : T_total)';
injection_gate = double(t_vector >= T_SETTLE);

for k = 1:3

    fprintf('---------------------------------------------------------\n');
    fprintf('Case %d: %s  (%.3f A)\n', k, labels{k}, amps(k));

    % ── INJECTION 1: q-sim = d-paper ─────────────────────────────────────
    assignin('base', 'id_perturb', ...
        timeseries(zeros(size(t_vector)), t_vector));
    assignin('base', 'iq_perturb', ...
        timeseries(amps(k)*sin(2*pi*fp_test*t_vector).*injection_gate, t_vector));

    simOut1 = sim(MODEL_NAME, [0 T_total], simOpts);
    [Vd1, Vq1, Id1, Iq1] = extract_phasors_local(simOut1, fp_test, Ts, T_SETTLE, N_samples);

    % Diagnostic plot first amplitude only
    if k == 1
        plot_diagnostic(simOut1, T_SETTLE, 10, 'Inj1 q-sim/d-paper');
        fprintf('  DC check: Id=%.3f A  Iq=%.3f A  Vd=%.3f V  Vq=%.3f V\n', ...
            mean(extract_column(simOut1.Id_Iq_data,1)), ...
            mean(extract_column(simOut1.Id_Iq_data,2)), ...
            mean(extract_column(simOut1.Vd_Vq_data,1)), ...
            mean(extract_column(simOut1.Vd_Vq_data,2)));
    end

    % ── INJECTION 2: d-sim = q-paper ─────────────────────────────────────
    assignin('base', 'id_perturb', ...
        timeseries(amps(k)*sin(2*pi*fp_test*t_vector).*injection_gate, t_vector));
    assignin('base', 'iq_perturb', ...
        timeseries(zeros(size(t_vector)), t_vector));

    simOut2 = sim(MODEL_NAME, [0 T_total], simOpts);
    [Vd2, Vq2, Id2, Iq2] = extract_phasors_local(simOut2, fp_test, Ts, T_SETTLE, N_samples);

    % ── Phasor dump ───────────────────────────────────────────────────────
    fprintf('  Phasors:\n');
    fprintf('    Iq1=%+.4f%+.4fj  Vq1=%+.4f%+.4fj  Vd1=%+.4f%+.4fj\n', ...
        real(Iq1),imag(Iq1), real(Vq1),imag(Vq1), real(Vd1),imag(Vd1));
    fprintf('    Id2=%+.4f%+.4fj  Vd2=%+.4f%+.4fj  Vq2=%+.4f%+.4fj\n', ...
        real(Id2),imag(Id2), real(Vd2),imag(Vd2), real(Vq2),imag(Vq2));

    % ── Matrix (paper frame: row1=d-paper=q-sim, row2=q-paper=d-sim) ─────
    %
    %  Paper eq.(19):  Y = [V_d1  V_d2]^-1 * [I_d1  I_d2]
    %                      [V_q1  V_q2]        [I_q1  I_q2]
    %
    %  In paper frame:  d-paper = q-sim,  q-paper = d-sim
    %  So:  V_d(paper) = Vq(sim),   I_d(paper) = Iq(sim)
    %       V_q(paper) = Vd(sim),   I_q(paper) = Id(sim)
    % ─────────────────────────────────────────────────────────────────────
    I_mat = [Iq1, Iq2;   % row 1: d-paper current  = q-sim current
             Id1, Id2];  % row 2: q-paper current  = d-sim current

    V_mat = [Vq1, Vq2;   % row 1: d-paper voltage  = q-sim voltage
             Vd1, Vd2];  % row 2: q-paper voltage  = d-sim voltage

    rc = rcond(V_mat);
    fprintf('  rcond(V_mat) = %.3e\n', rc);
    if rc < 1e-6
        warning('V_mat near-singular — skipping k=%d.', k);
        Z_results(:,k) = NaN;
        clear simOut1 simOut2;
        continue
    end

    Y_sim = V_mat \ I_mat;   % admittance, paper frame
    Z_sim = inv(Y_sim);      % impedance  — what Fig.9 plots

    Zdd = Z_sim(1,1);
    Zdq = Z_sim(1,2);
    Zqd = Z_sim(2,1);
    Zqq = Z_sim(2,2);

    fprintf('  Impedance (paper frame, dB):\n');
    fprintf('    Zdd: %.4f Ohm  ang=%+.1f deg  %+.1f dB\n', ...
        abs(Zdd), angle(Zdd)*180/pi, 20*log10(abs(Zdd)));
    fprintf('    Zdq: %.4f Ohm  ang=%+.1f deg  %+.1f dB\n', ...
        abs(Zdq), angle(Zdq)*180/pi, 20*log10(abs(Zdq)));
    fprintf('    Zqd: %.4f Ohm  ang=%+.1f deg  %+.1f dB\n', ...
        abs(Zqd), angle(Zqd)*180/pi, 20*log10(abs(Zqd)));
    fprintf('    Zqq: %.4f Ohm  ang=%+.1f deg  %+.1f dB\n', ...
        abs(Zqq), angle(Zqq)*180/pi, 20*log10(abs(Zqq)));

    Z_results(:,k) = [Zdd; Zdq; Zqd; Zqq];

    clear simOut1 simOut2;
end

% ── Pass/fail ─────────────────────────────────────────────────────────────
fprintf('\n=========================================================\n');
fprintf('Consistency (half vs double):\n');
diff_hd = abs(abs(Z_results(:,1)) - abs(Z_results(:,3))) ./ ...
          (abs(Z_results(:,2)) + eps) * 100;
names_Z = {'Zdd','Zdq','Zqd','Zqq'};
for i = 1:4
    fprintf('  %s: %.2f%%\n', names_Z{i}, diff_hd(i));
end
if max(diff_hd) < 5
    fprintf('PASS\n');
else
    fprintf('FAIL — max variation %.2f%%\n', max(diff_hd));
end

close_system(MODEL_NAME, 0);


% =========================================================================
% HELPERS
% =========================================================================
function plot_diagnostic(simOut, T_SETTLE, fig_num, ttl)
    Id_t = extract_column(simOut.Id_Iq_data, 1);
    Iq_t = extract_column(simOut.Id_Iq_data, 2);
    Vd_t = extract_column(simOut.Vd_Vq_data, 1);
    Vq_t = extract_column(simOut.Vd_Vq_data, 2);
    t    = simOut.tout;
    figure(fig_num); clf;
    subplot(2,2,1); plot(t,Id_t); xline(T_SETTLE,'r--'); title('Id'); grid on;
    subplot(2,2,2); plot(t,Iq_t); xline(T_SETTLE,'r--'); title('Iq'); grid on;
    subplot(2,2,3); plot(t,Vd_t); xline(T_SETTLE,'r--'); title('Vd'); grid on;
    subplot(2,2,4); plot(t,Vq_t); xline(T_SETTLE,'r--'); title('Vq'); grid on;
    sgtitle(ttl);
end

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

function [Vd_ph, Vq_ph, Id_ph, Iq_ph] = extract_phasors_local(simOut, fp, Ts, T_settle, N_samples)
    t      = simOut.tout;
    Id_all = extract_column(simOut.Id_Iq_data, 1);
    Iq_all = extract_column(simOut.Id_Iq_data, 2);
    Vd_all = extract_column(simOut.Vd_Vq_data, 1);
    Vq_all = extract_column(simOut.Vd_Vq_data, 2);

    i_start = find(t >= T_settle, 1, 'first');
    i_end   = i_start + N_samples - 1;
    if i_end > length(t)
        warning('Truncating extraction window.');
        i_end     = length(t);
        N_samples = i_end - i_start + 1;
    end

    w     = hann(N_samples);
    w_sum = sum(w);
    scale = 2 / w_sum;

    Vd_ph = fft((Vd_all(i_start:i_end)-mean(Vd_all(i_start:i_end))).*w, 2^nextpow2(N_samples*4));
    Vq_ph = fft((Vq_all(i_start:i_end)-mean(Vq_all(i_start:i_end))).*w, 2^nextpow2(N_samples*4));
    Id_ph = fft((Id_all(i_start:i_end)-mean(Id_all(i_start:i_end))).*w, 2^nextpow2(N_samples*4));
    Iq_ph = fft((Iq_all(i_start:i_end)-mean(Iq_all(i_start:i_end))).*w, 2^nextpow2(N_samples*4));

    NFFT  = 2^nextpow2(N_samples*4);
    fs    = 1/Ts;
    f_ax  = (0:NFFT-1)*fs/NFFT;
    [~,k_fp] = min(abs(f_ax - fp));

    Vd_ph = Vd_ph(k_fp) * scale;
    Vq_ph = Vq_ph(k_fp) * scale;
    Id_ph = Id_ph(k_fp) * scale;
    Iq_ph = Iq_ph(k_fp) * scale;
end