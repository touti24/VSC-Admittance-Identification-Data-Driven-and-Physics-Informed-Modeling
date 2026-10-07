function GPLL = pll_tf(Kp_pll, Ki_pll, Vd_ss)
% =========================================================================
% pll_tf.m
%
% Builds the linearised SRF-PLL transfer function as a MATLAB tf object.
%
% From Eq. (3) of Zhang et al. 2023:
%
%   GPLL(s) = (Kp_pll * s + Ki_pll)
%             ─────────────────────────────────────────
%             s^2 + (Kp_pll * Vd_ss) * s + Ki_pll * Vd_ss
%
% DERIVATION NOTE:
%   The PLL PI controller is: C_pll(s) = Kp_pll + Ki_pll/s
%   The plant (voltage to phase) is: 1/s (integrator × Vd_ss)
%   The closed-loop is: C(s)*Vd_ss / (s + C(s)*Vd_ss)
%   which gives Eq. 3 after simplification.
%
% INPUTS:
%   Kp_pll  — PLL proportional gain
%   Ki_pll  — PLL integral gain
%   Vd_ss   — steady-state d-axis voltage at PCC [V]
%              (= sqrt(2)*220 = 311 V for PLL-aligned operation)
%
% OUTPUT:
%   GPLL    — MATLAB tf object (Control System Toolbox)
%
% USAGE:
%   GPLL = pll_tf(1.2, 257, 311);
%   bode(GPLL);    % check PLL bandwidth
% =========================================================================

    num = [Kp_pll,  Ki_pll];
    den = [1,  Kp_pll * Vd_ss,  Ki_pll * Vd_ss];

    GPLL = tf(num, den);

    % Print PLL bandwidth for sanity check
    % Crossover frequency should be ~10-50 Hz for typical grid-connected VSC
    try
        bw = bandwidth(GPLL);
        fprintf('pll_tf: PLL bandwidth = %.1f Hz\n', bw / (2*pi));
    catch
        % bandwidth() may not be available in all toolbox versions
    end
end