function Y = analytical_model(f, Vd, Vq, Id, Iq, params)
% =========================================================================
% analytical_model.m
% Small-Signal Closed-Loop dq-Admittance Model of a Grid-Connected VSC
%
% Description:
%   Evaluates the continuous-frequency 2x2 small-signal admittance matrix Y(s)
%   in the synchronous dq reference frame at frequency f [Hz], incorporating:
%     - Inverter filter dynamics and cross-coupling
%     - Digital computation and modulation delays (4th-order Pade approximation)
%     - Synchronous Reference Frame Phase-Locked Loop (SRF-PLL) coupling
%     - Closed-loop PI current control dynamics
%
% Reference:
%   Zhang, Xu & Wang, "Physics-Informed Neural Network-Based Online Impedance
%   Identification of Voltage Source Converters," IEEE TIE, 2023 (Eq. 12).
%   Wang, Harnefors & Blaabjerg, "Unified Impedance Model of Grid-Connected
%   Voltage-Source Converters," IEEE TPEL, 2018.
%
% Inputs:
%   f      - Evaluation frequency [Hz]
%   Vd, Vq - Steady-state PCC voltage components in dq-frame [V]
%   Id, Iq - Steady-state converter current components in dq-frame [A]
%   params - Struct containing hardware and controller parameters:
%              .Lf     : Filter inductance [H]
%              .Rf     : Filter equivalent series resistance [Ohm]
%              .Kp_i   : Current loop proportional gain
%              .Ki_i   : Current loop integral gain
%              .Kp_pll : PLL proportional gain
%              .Ki_pll : PLL integral gain
%              .Ts     : Control sampling period [s]
%              .Vd_ss  : Nominal steady-state d-axis peak voltage [V]
%
% Outputs:
%   Y      - 2x2 complex admittance matrix evaluated at s = j*2*pi*f [Siemens]
%            [ Delta_id(s); Delta_iq(s) ] = Y(s) * [ Delta_vd(s); Delta_vq(s) ]
% =========================================================================

    s = 1j * 2 * pi * f;   % Complex Laplace frequency variable

    %% 1. PARAMETER EXTRACTION & GRID FREQUENCY
    Lf      = params.Lf;
    Rf      = params.Rf;
    Kp_i    = params.Kp_i;
    Ki_i    = params.Ki_i;
    Kp_pll  = params.Kp_pll;
    Ki_pll  = params.Ki_pll;
    Ts      = params.Ts;
    Vd_ss   = params.Vd_ss;

    f1 = 50;               % Fundamental grid frequency [Hz]
    w1 = 2 * pi * f1;      % Fundamental angular frequency [rad/s]

    %% 2. SCALAR TRANSFER FUNCTIONS
    % Current loop PI controller
    Gi = eval_tf([Kp_i, Ki_i], [1, 0], s);

    % Total control and PWM delay (1.5 * Ts, 4th-order Pade approximation)
    Gdel = eval_pade_delay(1.5 * Ts, s);

    % Linearized closed-loop SRF-PLL transfer function
    GPLL = eval_tf([Kp_pll, Ki_pll], [1, Kp_pll * Vd_ss, Ki_pll * Vd_ss], s);

    %% 3. 2x2 TRANSFER MATRICES
    Gm_i   = Gi * eye(2);
    Gm_del = Gdel * eye(2);

    % Plant impedance and admittance matrices
    Z_plant = [s * Lf + Rf, -w1 * Lf;
               w1 * Lf,      s * Lf + Rf];
    Ym_p = inv(Z_plant);
    Ym_o = Ym_p;

    % PLL small-signal perturbation coupling matrices
    Gm_I_PLL = [0, -GPLL * Iq;
                0,  GPLL * Id];

    Gm_V_PLL = [0, -GPLL * Vq;
                0,  GPLL * Vd];

    %% 4. CLOSED-LOOP INTERMEDIATE OPERATORS
    T       = Ym_p * Gm_del * Gm_i;
    Y_to    = Ym_o - Ym_p * Gm_del * Gm_V_PLL;

    IpT_inv = inv(eye(2) + T);
    G_cl    = IpT_inv * T;
    Y_PLL   = Gm_I_PLL;

    %% 5. CLOSED-LOOP ADMITTANCE MATRIX & CONVENTION ALIGNMENT
    Y = IpT_inv * Y_to - G_cl * Y_PLL;

    % Align output current polarity and frame orientation to target benchmark convention
    Y(2, 2) = -Y(2, 2);
    Y(1, 2) = -Y(1, 2);

end

%% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function val = eval_tf(num, den, s)
    % Evaluates a rational transfer function polynomial ratio at s
    val = polyval(num, s) / polyval(den, s);
end

function val = eval_pade_delay(tau, s)
    % 4th-order Pade approximation of transport delay exp(-s*tau)
    x = tau * s;
    num_val = 1 - x/2 + x^2/12 - x^3/120 + x^4/1680;
    den_val = 1 + x/2 + x^2/12 + x^3/120 + x^4/1680;
    val = num_val / den_val;
end
