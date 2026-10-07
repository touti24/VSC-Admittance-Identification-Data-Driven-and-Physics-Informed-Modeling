%% plot_bode_admittance.m
% =========================================================================
% Analytical dq-Admittance Bode Frequency Response Visualization (VSC I)
%
% Description:
%   Evaluates and renders the continuous small-signal frequency response
%   matrix Y(jw) for the baseline converter (VSC I) across the 1–100 Hz band.
%   Generates 2x4 Bode subplots (Magnitude in dB, Unwrapped Phase in degrees)
%   for all four admittance matrix elements:
%       [ Ydd(jw), Ydq(jw);
%         Yqd(jw), Yqq(jw) ]
%
% Reference:
%   Zhang, Xu & Wang, "Physics-Informed Neural Network-Based Online Impedance
%   Identification of Voltage Source Converters," IEEE TIE, 2023 (Eq. 12).
%   Hardware & control parameters from Table I (VSC I).
% =========================================================================

clear; clc; close all;

%% 1. HARDWARE & CONTROLLER SPECIFICATIONS (Table I)
params.Lf      = 1.0e-3;        % [H] Inverter-side filter inductance
params.Rf      = 0.3e-3;        % [Ohm] Equivalent series resistance
params.Kp_i    = 10.5;          % Current loop proportional gain
params.Ki_i    = 5900.0;        % Current loop integral gain [1/s]
params.Kp_pll  = 1.2;           % SRF-PLL proportional gain [rad/(V*s)]
params.Ki_pll  = 257.0;         % SRF-PLL integral gain [rad/(V*s^2)]
params.Ts      = 1 / 20000;     % [s] Control sampling period (fsw = 20 kHz)
params.Vd_ss   = 311.0;         % [V] Steady-state peak phase voltage

%% 2. OPERATING-POINT DEFINITION
Vd = 311.0;                     % [V] d-axis PCC peak voltage (PLL-aligned)
Vq = 0.0;                       % [V] q-axis PCC voltage
Id = 50.0;                      % [A] Active current setpoint
Iq = 0.0;                       % [A] Reactive current setpoint (unity PF)

%% 3. FREQUENCY RESPONSE EVALUATION (1–100 Hz)
f_vec = logspace(0, 2, 500);    % 500 logarithmic frequency points

% Pre-allocate frequency-response arrays
Ydd = zeros(1, length(f_vec));
Ydq = zeros(1, length(f_vec));
Yqd = zeros(1, length(f_vec));
Yqq = zeros(1, length(f_vec));

for k = 1:length(f_vec)
    Y = analytical_model(f_vec(k), Vd, Vq, Id, Iq, params);
    Ydd(k) = Y(1, 1);
    Ydq(k) = Y(1, 2);
    Yqd(k) = Y(2, 1);
    Yqq(k) = Y(2, 2);
end

%% 4. BODE DIAGRAM RENDERING
comps   = {Ydd, Ydq, Yqd, Yqq};
labels  = {'Y_{dd}', 'Y_{dq}', 'Y_{qd}', 'Y_{qq}'};
colors  = {'#0072BD', '#D95319', '#77AC30', '#7E2F8E'};

fig = figure('Name', 'dq-Admittance Bode - VSC I', ...
             'NumberTitle', 'off', ...
             'Position', [100, 80, 1100, 820]);

for idx = 1:4
    Y_comp  = comps{idx};
    mag_dB  = 20 * log10(abs(Y_comp));
    pha_deg = unwrap(angle(Y_comp)) * (180 / pi);

    % Magnitude Subplot (Top Row)
    ax_mag = subplot(2, 4, idx);
    semilogx(f_vec, mag_dB, 'Color', colors{idx}, 'LineWidth', 1.8);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(labels{idx}, 'FontWeight', 'bold');
    xlim([1, 100]);
    set(ax_mag, 'XTick', [1, 2, 5, 10, 20, 50, 100]);

    % Phase Subplot (Bottom Row)
    ax_pha = subplot(2, 4, idx + 4);
    semilogx(f_vec, pha_deg, 'Color', colors{idx}, 'LineWidth', 1.8);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Phase (deg)');
    title(labels{idx}, 'FontWeight', 'bold');
    xlim([1, 100]);
    set(ax_pha, 'XTick', [1, 2, 5, 10, 20, 50, 100]);
end

sgtitle(sprintf('dq-Admittance Frequency Response - VSC I (Analytical Model)\nV_d = %g V, V_q = %g V, I_d = %g A, I_q = %g A', ...
                Vd, Vq, Id, Iq), 'FontSize', 11, 'FontWeight', 'bold');

%% 5. FIGURE EXPORT
out_dir   = fileparts(mfilename('fullpath'));
save_name = sprintf('bode_admittance_VSC1_Id%dA', round(Id));
saveas(fig, fullfile(out_dir, [save_name, '.png']));
saveas(fig, fullfile(out_dir, [save_name, '.fig']));
fprintf('[plot_bode_admittance] Figures saved: %s.{png,fig}\n', save_name);

%% 6. SPOT-CHECK VERIFICATION TABLE
fprintf('\n--- Spot-Check Numerical Verification (Id = %g A) ---\n', Id);
for f_check = [1, 10, 50, 100]
    Y_chk = analytical_model(f_check, Vd, Vq, Id, Iq, params);
    fprintf('f = %3d Hz | Ydd: %+6.2f dB <%+7.2f deg | Yqq: %+6.2f dB <%+7.2f deg\n', ...
        f_check, ...
        20 * log10(abs(Y_chk(1, 1))), angle(Y_chk(1, 1)) * 180 / pi, ...
        20 * log10(abs(Y_chk(2, 2))), angle(Y_chk(2, 2)) * 180 / pi);
end
