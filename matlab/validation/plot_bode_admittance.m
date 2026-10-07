% =========================================================================
% plot_bode_admittance.m
%
% Plots Bode magnitude and phase of all four dq-admittance components
% (Ydd, Ydq, Yqd, Yqq) for VSC I parameters at a representative
% operating point (Id = 50 A, Iq = 0 A, unity power factor, PLL-aligned).
%
% Uses analytical_model.m which implements Eq. 12 of:
%   Zhang, Xu & Wang, IEEE Trans. Ind. Electron., Vol. 70 No. 4, 2023
%
% Operating-point choices (underdetermined in paper — see thesis notes):
%   Iq  = 0   : unity power factor (standard grid injection mode)
%   Vq  = 0   : PLL alignment sets q-axis voltage to zero in steady state
%   Vd  = 311 V: peak of 220 V_RMS line-to-neutral grid voltage
%   Id  = 50 A : mid-range representative value (paper sweeps 20-70 A)
%
% VSC I parameters from Table I of Zhang et al. 2023:
%   fsw = 20 kHz  →  Ts = 1/20000 s
%   Lf  = 1 mH,  Rf = 0.3 mΩ,  Kp_i = 10.5,  Ki_i = 5900
%   Kp_pll = 1.2,  Ki_pll = 257
% =========================================================================

clear; clc; close all;

%% ── 1. VSC I PARAMETERS (Table I, Zhang et al. 2023) ────────────────────
params.Lf      = 1e-3;          % [H]   inverter-side inductor
params.Rf      = 0.3e-3;        % [Ohm] parasitic resistance (0.3 mΩ)
params.Kp_i    = 10.5;          % current controller proportional gain
params.Ki_i    = 5900;          % current controller integral gain [1/s]
params.Kp_pll  = 1.2;           % PLL proportional gain [rad/s/V]
params.Ki_pll  = 257;           % PLL integral gain [rad/s²/V]
params.Ts      = 1/20000;       % [s]   sampling period (fsw = 20 kHz)
params.Vd_ss   = 311;           % [V]   steady-state d-axis voltage (220V_RMS peak)

%% ── 2. REPRESENTATIVE OPERATING POINT ───────────────────────────────────
Vd = 311;    % [V]  d-axis PCC voltage (PLL-aligned)
Vq = 0;      % [V]  q-axis PCC voltage (= 0 at PLL alignment)
Id = 50;     % [A]  d-axis current (mid-range, paper sweeps 20-70 A)
Iq = 0;      % [A]  q-axis current (unity power factor)

%% ── 3. FREQUENCY SWEEP ───────────────────────────────────────────────────
% Match the paper's measurement range: 1 to 100 Hz.
% Use fine logarithmic spacing for smooth Bode curves.
f_vec = logspace(0, 2, 500);   % 1 Hz to 100 Hz, 500 points

% Pre-allocate storage for all four Ykl components
Ydd = zeros(1, length(f_vec));
Ydq = zeros(1, length(f_vec));
Yqd = zeros(1, length(f_vec));
Yqq = zeros(1, length(f_vec));

for k = 1:length(f_vec)
    Y = analytical_model(f_vec(k), Vd, Vq, Id, Iq, params);
    Ydd(k) = Y(1,1);
    Ydq(k) = Y(1,2);
    Yqd(k) = Y(2,1);
    Yqq(k) = Y(2,2);
end

%% ── 4. BODE PLOTS ────────────────────────────────────────────────────────
% Convert to dB for magnitude; degrees for phase.
% Phase is unwrapped to avoid ±180° jumps in display.
comps   = {Ydd,  Ydq,  Yqd,  Yqq};
labels  = {'Y_{dd}', 'Y_{dq}', 'Y_{qd}', 'Y_{qq}'};
colors  = {'#0072BD', '#D95319', '#77AC30', '#7E2F8E'};

fig = figure('Name', 'dq-Admittance Bode — VSC I', ...
             'NumberTitle', 'off', ...
             'Position', [100 80 1100 820]);

for idx = 1:4
    Y_comp = comps{idx};
    mag_dB = 20 * log10(abs(Y_comp));
    pha_deg = unwrap(angle(Y_comp)) * (180/pi);

    % Magnitude subplot (top row)
    ax_mag = subplot(2, 4, idx);
    semilogx(f_vec, mag_dB, 'Color', colors{idx}, 'LineWidth', 1.8);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(labels{idx});
    xlim([1 100]);
    set(ax_mag, 'XTick', [1 2 5 10 20 50 100]);

    % Phase subplot (bottom row)
    ax_pha = subplot(2, 4, idx+4);
    semilogx(f_vec, pha_deg, 'Color', colors{idx}, 'LineWidth', 1.8);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Phase (deg)');
    title(labels{idx});
    xlim([1 100]);
    set(ax_pha, 'XTick', [1 2 5 10 20 50 100]);
end

% Super-title with operating point annotation
sgtitle(sprintf(['dq-Admittance Bode — VSC I (Analytical, Eq. 12)\n' ...
                 'V_d = %g V, V_q = %g V, I_d = %g A, I_q = %g A'], ...
                 Vd, Vq, Id, Iq), 'FontSize', 11, 'FontWeight', 'bold');

%% ── 5. SAVE FIGURE ───────────────────────────────────────────────────────
% Save as both PNG and .fig for inclusion in report / further editing.
save_name = sprintf('bode_admittance_VSC1_Id%dA', Id);
saveas(fig, fullfile(fileparts(mfilename('fullpath')), [save_name '.png']));
saveas(fig, fullfile(fileparts(mfilename('fullpath')), [save_name '.fig']));
fprintf('[plot_bode_admittance] Figures saved: %s.{png,fig}\n', save_name);

%% ── 6. PRINT SPOT VALUES FOR SANITY CHECK ───────────────────────────────
% At f=1 Hz the admittance should be dominated by the L-filter:
%   Ydd ≈ 1/(Rf + j*2pi*1*Lf) ≈ 1/0.3e-3 = 3333 S  (huge, Rf tiny)
% At f=100 Hz: |Y| ≈ 1/(2pi*100*1e-3) ≈ 1.59 S → -3.97 dB  (L dominates)
fprintf('\n--- Spot-check values (VSC I, Id=%dA) ---\n', Id);
for f_check = [1, 10, 50, 100]
    Y_chk = analytical_model(f_check, Vd, Vq, Id, Iq, params);
    fprintf('f=%4g Hz | Ydd: %+6.2f dB ∠%+7.2f°  | Yqq: %+6.2f dB ∠%+7.2f°\n', ...
        f_check, ...
        20*log10(abs(Y_chk(1,1))), angle(Y_chk(1,1))*180/pi, ...
        20*log10(abs(Y_chk(2,2))), angle(Y_chk(2,2))*180/pi);
end
