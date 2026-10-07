function GPLL = pll_tf(Kp_pll, Ki_pll, Vd_ss)
% =========================================================================
% pll_tf.m
% Linearized Synchronous Reference Frame PLL (SRF-PLL) Transfer Function
%
% Description:
%   Constructs the closed-loop small-signal transfer function of an SRF-PLL:
%
%                 Kp_pll * s + Ki_pll
%     GPLL(s) = --------------------------------------------
%               s^2 + (Kp_pll * Vd_ss) * s + (Ki_pll * Vd_ss)
%
% Reference:
%   Zhang, Xu & Wang, "Physics-Informed Neural Network-Based Online Impedance
%   Identification of Voltage Source Converters," IEEE TIE, 2023 (Eq. 3).
%
% Inputs:
%   Kp_pll - Proportional gain of the PLL PI controller [rad/(V*s)]
%   Ki_pll - Integral gain of the PLL PI controller [rad/(V*s^2)]
%   Vd_ss  - Steady-state d-axis voltage at PCC [V] (nominal peak ~ 311 V)
%
% Outputs:
%   GPLL   - Linear continuous-time transfer function (MATLAB tf object)
% =========================================================================

    %% 1. Polynomial Formulation
    num = [Kp_pll, Ki_pll];
    den = [1, Kp_pll * Vd_ss, Ki_pll * Vd_ss];

    %% 2. Transfer Function Construction
    GPLL = tf(num, den);

    %% 3. Bandwidth Verification
    try
        bw = bandwidth(GPLL);
        fprintf('[pll_tf] Closed-loop PLL bandwidth: %.2f Hz\n', bw / (2 * pi));
    catch
        % Fallback if Control System Toolbox bandwidth() is unavailable
    end

end
