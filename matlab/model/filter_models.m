function [Yp, Yo] = filter_models(filter_type, params)
% =========================================================================
% filter_models.m
% Continuous-Time Converter and Grid Interface Filter Transfer Functions
%
% Description:
%   Constructs continuous-time admittance transfer functions Yp(s) and Yo(s)
%   for small-signal impedance/admittance modeling of grid-connected VSCs:
%     - 'L'   : Single-inductor interface filter (standard L-filter).
%     - 'LCL' : Third-order LCL filter with capacitor current feedback
%               active damping.
%
% Reference:
%   Zhang, Xu & Wang, "Physics-Informed Neural Network-Based Online Impedance
%   Identification of Voltage Source Converters," IEEE TIE, 2023.
%   Pan et al., "An Improved Capacitor-Current Active Damping for LCL-Filter,"
%   IEEE TPEL, 2015.
%
% Inputs:
%   filter_type - String specifying filter topology: 'L' or 'LCL'
%   params      - Struct containing hardware parameters:
%                   For 'L':
%                     .Lf : Inverter-side filter inductance [H]
%                     .Rf : Filter equivalent series resistance [Ohm]
%                   For 'LCL':
%                     .Lf : Inverter-side filter inductance [H]
%                     .Rf : Inverter-side series resistance [Ohm]
%                     .L2 : Grid-side filter inductance [H]
%                     .Cf : Filter capacitance [F]
%                     .Kd : Active damping feedback gain [Ohm]
%
% Outputs:
%   Yp          - Inverter-side admittance transfer function Yp(s) (MATLAB tf)
%   Yo          - Grid-side (PCC) admittance transfer function Yo(s) (MATLAB tf)
% =========================================================================

    switch lower(string(filter_type))

        %% 1. L-FILTER TOPOLOGY
        case "l"
            Lf = params.Lf;
            Rf = params.Rf;

            % Converter-side and grid-side admittances are identical
            % Yp(s) = Yo(s) = 1 / (Lf*s + Rf)
            Yp = tf(1, [Lf, Rf]);
            Yo = Yp;

            fprintf("[filter_models] L-Filter initialized: Lf = %.2f mH, Rf = %.3f mOhm\n", ...
                    Lf * 1e3, Rf * 1e3);

        %% 2. LCL-FILTER TOPOLOGY (WITH ACTIVE DAMPING)
        case "lcl"
            Lf = params.Lf;     % Inverter-side inductance [H]
            Rf = params.Rf;     % Inverter-side series resistance [Ohm]
            L2 = params.L2;     % Grid-side inductance [H]
            Cf = params.Cf;     % Shunt capacitance [F]
            Kd = params.Kd;     % Capacitor-current active damping gain [Ohm]

            % Inverter-side admittance with active damping emulation
            % Yp(s) = 1 / (Lf*s + Rf + Kd)
            Yp = tf(1, [Lf, Rf + Kd]);

            % Grid-side admittance (un-damped inductive grid interface)
            % Yo(s) = 1 / (L2*s)
            Yo = tf(1, [L2, 0]);

            fprintf("[filter_models] LCL-Filter initialized: Lf = %.2f mH, L2 = %.2f mH, Cf = %.1f uF, Kd = %.1f Ohm\n", ...
                    Lf * 1e3, L2 * 1e3, Cf * 1e6, Kd);

        otherwise
            error("[filter_models] Unknown filter topology '%s'. Valid options are 'L' or 'LCL'.", filter_type);
    end

end
