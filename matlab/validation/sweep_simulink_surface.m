% =========================================================================
% sweep_admittance_VSC_I.m
%
%  Reproduces Fig. 9 of the paper:
%    "Measured admittance dataset of VSC I"
%
%  Sweeps frequency fp = 1:1:100 Hz and d-axis current Id = 20:1:70 A
%  At each (fp, Id) point:
%    - Two injections (q-sim/d-paper, then d-sim/q-paper)
%    - Extract phasors via Hann-windowed FFT
%    - Build Y matrix in paper frame
%  Plots 8 subplots (magnitude + phase) for Ydd, Ydq, Yqd, Yqq
%
%  Convention:
%    paper d-axis = simulation q-axis
%    paper q-axis = simulation d-axis
%    Fig.9 plots impedance Y (not admittance Y) in dB
%
%  CHECKPOINT / RESUME
%    After every Id row the partial results are saved to
%    sweep_checkpoint_VSC_I.mat.  If the script is interrupted and
%    re-run it will automatically resume from the last completed row.
% =========================================================================

CHECKPOINT_FILE = 'sweep_checkpoint_VSC_I.mat';
FINAL_FILE      = 'sweep_results_VSC_I.mat';

MODEL_NAME = 'VSC_Iavrg';
run('params_VSC_I.m');       % loads Ts and VSC parameters

% ── Sweep axes (match paper Fig.9) ───────────────────────────────────────
fp_vec = 1 : 1 : 100;        % [Hz]   frequency axis
Id_vec = 20 : 1 : 70;        % [A]    d-axis current axis
Iq_ref = 0;                  % [A]    q-axis always zero

N_fp = length(fp_vec);
N_Id = length(Id_vec);

% ── Checkpoint: resume or start fresh ────────────────────────────────────
if isfile(CHECKPOINT_FILE)
    fprintf('>>> Checkpoint found — resuming from last saved row.\n');
    load(CHECKPOINT_FILE, ...
        'Ydd_mag','Ydd_ang','Ydq_mag','Ydq_ang', ...
        'Yqd_mag','Yqd_ang','Yqq_mag','Yqq_ang', ...
        'i_Id_done');
    i_Id_start = i_Id_done + 1;
    fprintf('    Resuming at Id index %d  (Id = %d A)\n\n', ...
        i_Id_start, Id_vec(i_Id_start));
else
    % Pre-allocate results  (rows=Id, cols=fp)
    Ydd_mag = zeros(N_Id, N_fp);   Ydd_ang = zeros(N_Id, N_fp);
    Ydq_mag = zeros(N_Id, N_fp);   Ydq_ang = zeros(N_Id, N_fp);
    Yqd_mag = zeros(N_Id, N_fp);   Yqd_ang = zeros(N_Id, N_fp);
    Yqq_mag = zeros(N_Id, N_fp);   Yqq_ang = zeros(N_Id, N_fp);
    i_Id_start = 1;
end

% ── Simulation settings ───────────────────────────────────────────────────
T_SETTLE   = 8.0;            % [s] settling time (same as sanity check)
N_cycles   = 5;              % perturbation cycles per injection
AMP_FRAC   = 0.05;           % 5% of operating point

simOpts = simset('SrcWorkspace', 'base');
load_system(MODEL_NAME);

total_sims = 2 * N_fp * N_Id;
sim_count  = 2 * N_fp * (i_Id_start - 1);   % already done
t_start    = tic;

fprintf('=========================================================\n');
fprintf('VSC I admittance sweep\n');
fprintf('  fp : %d to %d Hz  (%d points)\n', fp_vec(1),   fp_vec(end),   N_fp);
fprintf('  Id : %d to %d A   (%d points)\n', Id_vec(1),   Id_vec(end),   N_Id);
fprintf('  Total simulations: %d\n', total_sims);
if i_Id_start > 1
    fprintf('  Skipping first %d Id rows (already in checkpoint).\n', i_Id_start-1);
end
fprintf('=========================================================\n');

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
        assignin('base', 'id_perturb', ...
            timeseries(zeros(size(t_vector)), t_vector));
        assignin('base', 'iq_perturb', ...
            timeseries(amp * sin(2*pi*fp*t_vector) .* gate, t_vector));

        simOut1 = sim(MODEL_NAME, [0 T_total], simOpts);
        [Vd1,Vq1,Id1,Iq1] = extract_phasors_local(simOut1, fp, Ts, T_SETTLE, N_samples);
        clear simOut1;

        sim_count = sim_count + 1;

        % ── Injection 2: d-sim = q-paper ─────────────────────────────────
        assignin('base', 'id_perturb', ...
            timeseries(amp * sin(2*pi*fp*t_vector) .* gate, t_vector));
        assignin('base', 'iq_perturb', ...
            timeseries(zeros(size(t_vector)), t_vector));

        simOut2 = sim(MODEL_NAME, [0 T_total], simOpts);
        [Vd2,Vq2,Id2,Iq2] = extract_phasors_local(simOut2, fp, Ts, T_SETTLE, N_samples);
        clear simOut2;

        sim_count = sim_count + 1;

        % ── Build matrices (paper frame) ──────────────────────────────────
        %  row1 = d-paper = q-sim
        %  row2 = q-paper = d-sim
        I_mat = [Iq1, Iq2;
                 Id1, Id2];
        V_mat = [Vq1, Vq2;
                 Vd1, Vd2];

        rc = rcond(V_mat);
        if rc < 1e-10
            warning('Near-singular at fp=%.0f Hz, Id=%.0f A — filling NaN.', fp, Id_test);
            Ydd_mag(i_Id,i_fp) = NaN;  Ydd_ang(i_Id,i_fp) = NaN;
            Ydq_mag(i_Id,i_fp) = NaN;  Ydq_ang(i_Id,i_fp) = NaN;
            Yqd_mag(i_Id,i_fp) = NaN;  Yqd_ang(i_Id,i_fp) = NaN;
            Yqq_mag(i_Id,i_fp) = NaN;  Yqq_ang(i_Id,i_fp) = NaN;
            continue
        end

        Y_sim = V_mat \ I_mat;
        Y_mat = inv(Y_sim);

        Ydd_mag(i_Id,i_fp) = abs(Y_mat(1,1));
        Ydd_ang(i_Id,i_fp) = angle(Y_mat(1,1)) * 180/pi;
        Ydq_mag(i_Id,i_fp) = abs(Y_mat(1,2));
        Ydq_ang(i_Id,i_fp) = angle(Y_mat(1,2)) * 180/pi;
        Yqd_mag(i_Id,i_fp) = abs(Y_mat(2,1));
        Yqd_ang(i_Id,i_fp) = angle(Y_mat(2,1)) * 180/pi;
        Yqq_mag(i_Id,i_fp) = abs(Y_mat(2,2));
        Yqq_ang(i_Id,i_fp) = angle(Y_mat(2,2)) * 180/pi;

    end % fp loop

    % ── Checkpoint save after every Id row ───────────────────────────────
    i_Id_done = i_Id;
    save(CHECKPOINT_FILE, ...
        'Ydd_mag','Ydd_ang','Ydq_mag','Ydq_ang', ...
        'Yqd_mag','Yqd_ang','Yqq_mag','Yqq_ang', ...
        'i_Id_done');

    % Progress report after each Id row
    elapsed = toc(t_start);
    done_frac = (i_Id - i_Id_start + 1) / (N_Id - i_Id_start + 1);
    eta = elapsed / done_frac * (1 - done_frac);
    fprintf('  Id=%2d A done  (%d/%d)  elapsed=%s  ETA=%s  [checkpoint saved]\n', ...
        Id_test, i_Id, N_Id, ...
        format_time(elapsed), format_time(eta));

end % Id loop

close_system(MODEL_NAME, 0);

fprintf('\nSweep complete in %s\n', format_time(toc(t_start)));

% ── Save final workspace & remove checkpoint ─────────────────────────────
save(FINAL_FILE, ...
    'fp_vec','Id_vec', ...
    'Ydd_mag','Ydd_ang','Ydq_mag','Ydq_ang', ...
    'Yqd_mag','Yqd_ang','Yqq_mag','Yqq_ang');
fprintf('Results saved to %s\n\n', FINAL_FILE);

if isfile(CHECKPOINT_FILE)
    delete(CHECKPOINT_FILE);
    fprintf('Checkpoint file removed.\n\n');
end

% ── Plot Fig. 9 ───────────────────────────────────────────────────────────
plot_fig9(fp_vec, Id_vec, ...
    Ydd_mag, Ydd_ang, Ydq_mag, Ydq_ang, ...
    Yqd_mag, Yqd_ang, Yqq_mag, Yqq_ang);


% =========================================================================
% PLOT  —  reproduces Fig. 9 layout exactly
%   Top row:    magnitude in dB  (Ydd | Ydq | Yqd | Yqq)
%   Bottom row: phase in degrees (Ydd | Ydq | Yqd | Yqq)
% =========================================================================
function plot_fig9(fp_vec, Id_vec, ...
        Ydd_mag, Ydd_ang, Ydq_mag, Ydq_ang, ...
        Yqd_mag, Yqd_ang, Yqq_mag, Yqq_ang)

    [FP, ID] = meshgrid(fp_vec, Id_vec);

    labels  = {'Y_{dd}', 'Y_{dq}', 'Y_{qd}', 'Y_{qq}'};
    mag_all = {20*log10(Ydd_mag), 20*log10(Ydq_mag), ...
               20*log10(Yqd_mag), 20*log10(Yqq_mag)};
    ang_all = {Ydd_ang, Ydq_ang, Yqd_ang, Yqq_ang};

    figure('Name','Fig. 9 — Measured admittance dataset of VSC I', ...
           'Position', [50 50 1400 700]);

    for col = 1:4

        % ── Magnitude (top row) ───────────────────────────────────────────
        subplot(2, 4, col);
        surf(FP, ID, mag_all{col}, 'EdgeColor','none');
        view([-37.5, 30]);
        xlabel('Frequency/(Hz)'); ylabel('Current/(A)'); zlabel('Magnitude/(dB)');
        title(labels{col});
        colormap(gca, parula);
        shading interp;
        grid on;
        set(gca,'XDir','reverse');   % paper has frequency axis reversed

        % ── Phase (bottom row) ────────────────────────────────────────────
        subplot(2, 4, col + 4);
        surf(FP, ID, ang_all{col}, 'EdgeColor','none');
        view([-37.5, 30]);
        xlabel('Frequency/(Hz)'); ylabel('Current/(A)'); zlabel('Phase/(deg)');
        title(labels{col});
        colormap(gca, parula);
        shading interp;
        grid on;
        set(gca,'XDir','reverse');

    end

    sgtitle('Measured admittance dataset of VSC I  (Fig. 9 reproduction)');
end


% =========================================================================
% PHASOR EXTRACTION
%   Hann window + zero-padded FFT → complex phasor at fp
%   scale = 2/sum(w) recovers true peak amplitude
% =========================================================================
function [Vd_ph, Vq_ph, Id_ph, Iq_ph] = extract_phasors_local(simOut, fp, Ts, T_settle, N_samples)

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
    w_sum = sum(w);
    scale = 2 / w_sum;
    NFFT  = 2^nextpow2(N_samples * 4);
    fs    = 1 / Ts;
    f_ax  = (0 : NFFT-1) * fs / NFFT;
    [~, k_fp] = min(abs(f_ax - fp));

    Vd_ph = fft((Vd_all(i_start:i_end)-mean(Vd_all(i_start:i_end))).*w, NFFT);
    Vq_ph = fft((Vq_all(i_start:i_end)-mean(Vq_all(i_start:i_end))).*w, NFFT);
    Id_ph = fft((Id_all(i_start:i_end)-mean(Id_all(i_start:i_end))).*w, NFFT);
    Iq_ph = fft((Iq_all(i_start:i_end)-mean(Iq_all(i_start:i_end))).*w, NFFT);

    Vd_ph = Vd_ph(k_fp) * scale;
    Vq_ph = Vq_ph(k_fp) * scale;
    Id_ph = Id_ph(k_fp) * scale;
    Iq_ph = Iq_ph(k_fp) * scale;
end


% =========================================================================
% HELPER: extract one column regardless of To-Workspace format
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


% =========================================================================
% HELPER: format seconds as h m s string
% =========================================================================
function s = format_time(sec)
    sec = round(sec);
    h   = floor(sec / 3600);
    m   = floor(mod(sec, 3600) / 60);
    s_  = mod(sec, 60);
    if h > 0
        s = sprintf('%dh %dm %ds', h, m, s_);
    elseif m > 0
        s = sprintf('%dm %ds', m, s_);
    else
        s = sprintf('%ds', s_);
    end
end