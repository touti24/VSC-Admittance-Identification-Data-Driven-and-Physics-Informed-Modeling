%% sweep_admittance_VSC_I.m
% =========================================================================
% Full-Grid Simulation Admittance Sweep for VSC I (Fig. 9 Reproduction)
%
% Description:
%   Performs an automated small-signal perturbation sweep across frequency
%   (fp = 1:1:100 Hz) and active operating current (Id = 20:1:70 A) using the
%   average-value model (VSC_Iavrg) in Simulink.
%
% Protocol:
%   - For each operating point (fp, Id):
%       1. Two orthogonal current injections (q-sim/d-paper, d-sim/q-paper).
%       2. Hann-windowed FFT phasor extraction of fundamental voltages/currents.
%       3. Multi-channel impedance/admittance matrix reconstruction.
%   - Plots 2x4 3D response surfaces (Magnitude [dB] and Phase [deg]) matching
%     Fig. 9 of Zhang et al., IEEE TIE (2023).
%
% Checkpoint & Resume:
%   Partial progress is automatically saved to 'sweep_checkpoint_VSC_I.mat'
%   after every Id row, allowing seamless resumption if interrupted.
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & SIMULATION SETUP
CHECKPOINT_FILE = 'sweep_checkpoint_VSC_I.mat';
FINAL_FILE      = 'sweep_results_VSC_I.mat';
MODEL_NAME      = 'VSC_Iavrg';

run('params_VSC_I.m');

% Sweep operational coordinates (matching Fig. 9 specifications)
fp_vec = 1 : 1 : 100;                 % [Hz] Frequency range (100 points)
Id_vec = 20 : 1 : 70;                 % [A] Active current range (51 points)
Iq_ref = 0;                           % [A] Reactive current setpoint (unity PF)

N_fp = length(fp_vec);
N_Id = length(Id_vec);

% Simulation timing and perturbation parameters
T_SETTLE   = 8.0;                     % Pre-perturbation settling time [s]
N_cycles   = 5;                       % FFT observation cycles per frequency
AMP_FRAC   = 0.05;                    % 5% nominal perturbation amplitude

simOpts    = simset('SrcWorkspace', 'base');
load_system(MODEL_NAME);

%% 2. CHECKPOINT INITIALIZATION / RESUME LOGIC
if isfile(CHECKPOINT_FILE)
    fprintf('[Resume] Checkpoint detected. Resuming from last completed state...\n');
    load(CHECKPOINT_FILE, ...
        'Ydd_mag', 'Ydd_ang', 'Ydq_mag', 'Ydq_ang', ...
        'Yqd_mag', 'Yqd_ang', 'Yqq_mag', 'Yqq_ang', ...
        'i_Id_done');
    i_Id_start = i_Id_done + 1;
    fprintf('  • Resuming at Id index %d (Id = %d A)\n\n', i_Id_start, Id_vec(i_Id_start));
else
    % Pre-allocate response grids (Rows = Id, Columns = fp)
    Ydd_mag = zeros(N_Id, N_fp);   Ydd_ang = zeros(N_Id, N_fp);
    Ydq_mag = zeros(N_Id, N_fp);   Ydq_ang = zeros(N_Id, N_fp);
    Yqd_mag = zeros(N_Id, N_fp);   Yqd_ang = zeros(N_Id, N_fp);
    Yqq_mag = zeros(N_Id, N_fp);   Yqq_ang = zeros(N_Id, N_fp);
    i_Id_start = 1;
end

total_sims = 2 * N_fp * N_Id;
sim_count  = 2 * N_fp * (i_Id_start - 1);
t_start    = tic;

fprintf('=========================================================\n');
fprintf('VSC I Multi-Operating-Point Admittance Sweep\n');
fprintf('  Frequency Span : %d to %d Hz (%d points)\n', fp_vec(1), fp_vec(end), N_fp);
fprintf('  Current Span   : %d to %d A  (%d points)\n', Id_vec(1), Id_vec(end), N_Id);
fprintf('  Total Simulations: %d\n', total_sims);
fprintf('=========================================================\n\n');

%% 3. AUTOMATED GRID SWEEP EXECUTION
for i_Id = i_Id_start : N_Id
    Id_test = Id_vec(i_Id);
    amp     = AMP_FRAC * Id_test;

    assignin('base', 'Id_ref', Id_test);
    assignin('base', 'Iq_ref', Iq_ref);

    for i_fp = 1 : N_fp
        fp = fp_vec(i_fp);

        T_perturb = N_cycles / fp;
        N_samples = round(T_perturb / Ts);
        T_total   = T_SETTLE + T_perturb;
        t_vector  = (0 : Ts : T_total)';
        gate      = double(t_vector >= T_SETTLE);

        % ── Injection 1: q-sim = d-paper ─────────────────────────────────
        assignin('base', 'id_perturb', timeseries(zeros(size(t_vector)), t_vector));
        assignin('base', 'iq_perturb', timeseries(amp * sin(2*pi*fp*t_vector) .* gate, t_vector));

        simOut1 = sim(MODEL_NAME, [0, T_total], simOpts);
        [Vd1, Vq1, Id1, Iq1] = extract_phasors_local(simOut1, fp, Ts, T_SETTLE, N_samples);
        clear simOut1;

        sim_count = sim_count + 1;

        % ── Injection 2: d-sim = q-paper ─────────────────────────────────
        assignin('base', 'id_perturb', timeseries(amp * sin(2*pi*fp*t_vector) .* gate, t_vector));
        assignin('base', 'iq_perturb', timeseries(zeros(size(t_vector)), t_vector));

        simOut2 = sim(MODEL_NAME, [0, T_total], simOpts);
        [Vd2, Vq2, Id2, Iq2] = extract_phasors_local(simOut2, fp, Ts, T_SETTLE, N_samples);
        clear simOut2;

        sim_count = sim_count + 1;

        % ── Admittance Matrix Formulation ────────────────────────────────
        % Target frame: Row 1 = d-paper (q-sim), Row 2 = q-paper (d-sim)
        I_mat = [Iq1, Iq2;
                 Id1, Id2];
        V_mat = [Vq1, Vq2;
                 Vd1, Vd2];

        rc = rcond(V_mat);
        if rc < 1e-10
            warning('Ill-conditioned matrix at fp = %.1f Hz, Id = %.1f A. Writing NaN.', fp, Id_test);
            Ydd_mag(i_Id, i_fp) = NaN;  Ydd_ang(i_Id, i_fp) = NaN;
            Ydq_mag(i_Id, i_fp) = NaN;  Ydq_ang(i_Id, i_fp) = NaN;
            Yqd_mag(i_Id, i_fp) = NaN;  Yqd_ang(i_Id, i_fp) = NaN;
            Yqq_mag(i_Id, i_fp) = NaN;  Yqq_ang(i_Id, i_fp) = NaN;
            continue;
        end

        Y_sim = V_mat \ I_mat;
        Y_mat = inv(Y_sim);

        Ydd_mag(i_Id, i_fp) = abs(Y_mat(1, 1));
        Ydd_ang(i_Id, i_fp) = angle(Y_mat(1, 1)) * 180 / pi;
        Ydq_mag(i_Id, i_fp) = abs(Y_mat(1, 2));
        Ydq_ang(i_Id, i_fp) = angle(Y_mat(1, 2)) * 180 / pi;
        Yqd_mag(i_Id, i_fp) = abs(Y_mat(2, 1));
        Yqd_ang(i_Id, i_fp) = angle(Y_mat(2, 1)) * 180 / pi;
        Yqq_mag(i_Id, i_fp) = abs(Y_mat(2, 2));
        Yqq_ang(i_Id, i_fp) = angle(Y_mat(2, 2)) * 180 / pi;
    end

    % ── Save Incremental Row Checkpoint ──────────────────────────────────
    i_Id_done = i_Id;
    save(CHECKPOINT_FILE, ...
        'Ydd_mag', 'Ydd_ang', 'Ydq_mag', 'Ydq_ang', ...
        'Yqd_mag', 'Yqd_ang', 'Yqq_mag', 'Yqq_ang', ...
        'i_Id_done');

    elapsed   = toc(t_start);
    done_frac = (i_Id - i_Id_start + 1) / (N_Id - i_Id_start + 1);
    eta       = (elapsed / done_frac) * (1 - done_frac);
    fprintf('  • Id = %2d A completed (%2d/%2d) | Elapsed: %s | ETA: %s [Checkpoint saved]\n', ...
        Id_test, i_Id, N_Id, format_time(elapsed), format_time(eta));
end

close_system(MODEL_NAME, 0);

%% 4. FINAL WORKSPACE SERIALIZATION & CLEANUP
save(FINAL_FILE, ...
    'fp_vec', 'Id_vec', ...
    'Ydd_mag', 'Ydd_ang', 'Ydq_mag', 'Ydq_ang', ...
    'Yqd_mag', 'Yqd_ang', 'Yqq_mag', 'Yqq_ang');
fprintf('\n[Success] Final results saved -> %s\n', FINAL_FILE);

if isfile(CHECKPOINT_FILE)
    delete(CHECKPOINT_FILE);
    fprintf('[Cleanup] Temporary checkpoint file removed.\n');
end

%% 5. VISUALIZATION (FIGURE 9 REPRODUCTION)
plot_fig9(fp_vec, Id_vec, ...
    Ydd_mag, Ydd_ang, Ydq_mag, Ydq_ang, ...
    Yqd_mag, Yqd_ang, Yqq_mag, Yqq_ang);

%% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function plot_fig9(fp_vec, Id_vec, ...
        Ydd_mag, Ydd_ang, Ydq_mag, Ydq_ang, ...
        Yqd_mag, Yqd_ang, Yqq_mag, Yqq_ang)
    % Renders Fig. 9 response surfaces (Magnitude & Phase)
    [FP, ID] = meshgrid(fp_vec, Id_vec);

    labels  = {'Y_{dd}', 'Y_{dq}', 'Y_{qd}', 'Y_{qq}'};
    mag_all = {20 * log10(Ydd_mag), 20 * log10(Ydq_mag), ...
               20 * log10(Yqd_mag), 20 * log10(Yqq_mag)};
    ang_all = {Ydd_ang, Ydq_ang, Yqd_ang, Yqq_ang};

    figure('Name', 'Fig. 9 Reproduction - Measured Admittance of VSC I', ...
           'Position', [50, 50, 1400, 700]);

    for col = 1:4
        % Magnitude Subplot (Top Row)
        subplot(2, 4, col);
        surf(FP, ID, mag_all{col}, 'EdgeColor', 'none');
        view([-37.5, 30]);
        xlabel('Frequency/(Hz)'); ylabel('Current/(A)'); zlabel('Magnitude/(dB)');
        title(labels{col}, 'FontWeight', 'bold');
        colormap(gca, parula);
        shading interp;
        grid on;
        set(gca, 'XDir', 'reverse');

        % Phase Subplot (Bottom Row)
        subplot(2, 4, col + 4);
        surf(FP, ID, ang_all{col}, 'EdgeColor', 'none');
        view([-37.5, 30]);
        xlabel('Frequency/(Hz)'); ylabel('Current/(A)'); zlabel('Phase/(deg)');
        title(labels{col}, 'FontWeight', 'bold');
        colormap(gca, parula);
        shading interp;
        grid on;
        set(gca, 'XDir', 'reverse');
    end

    sgtitle('Measured Admittance Dataset of VSC I (Fig. 9 Reproduction)', 'FontWeight', 'bold');
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

function s = format_time(sec)
    % Formats duration in seconds into human-readable string
    sec = round(sec);
    h   = floor(sec / 3600);
    m   = floor(mod(sec, 3600) / 60);
    s_  = mod(sec, 60);
    if h > 0
        s = sprintf('%dh %02dm %02ds', h, m, s_);
    elseif m > 0
        s = sprintf('%dm %02ds', m, s_);
    else
        s = sprintf('%ds', s_);
    end
end
