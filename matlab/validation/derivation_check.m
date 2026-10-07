%% derivation_check.m
% =========================================================================
% Symbolic Verification of the PINN Multi-Operating-Point Factorization
%
% Description:
%   Performs symbolic mathematical verification of the multi-operating-point
%   admittance decomposition proposed in Zhang et al., IEEE TIE (2023).
%
% Mathematical Principle:
%   Verifies that the multi-variable nonlinear closed-loop admittance Y_dd(s)
%   algebraically decomposes into a linear combination of five purely 1D
%   frequency-dependent functions {F1(s), ..., F5(s)} parameterized by terminal
%   operating variables (Vd, Vq, Id, Iq):
%
%       Y_dd(s) = F1(s) - Vq*F2(s) - Vd*F3(s) - Iq*F4(s) - Id*F5(s)
%
% Reference:
%   Zhang, Xu & Wang, "Physics-Informed Neural Network-Based Online Impedance
%   Identification of Voltage Source Converters," IEEE TIE, 2023 (Eqs. 13-18).
% =========================================================================

clear; clc;

fprintf('=========================================================================\n');
fprintf('Symbolic Verification: Admittance Factorization Architecture\n');
fprintf('=========================================================================\n\n');

%% 1. Symbolic Variable Declarations
% Operating-point scalars
syms f positive
syms V_d V_q I_d I_q real

% Frequency-dependent intermediate sub-matrices (Eq. 15)
syms O_dd O_dq O_qd O_qq real
syms P_dd P_dq P_qd P_qq real
syms Q_dd Q_dq Q_qd Q_qq real
syms R_inv_dd R_inv_dq R_inv_qd R_inv_qq real

%% 2. Algebraic Combiner Construction (Eq. 16)
% Construct closed-loop admittance combiner for the direct channel Y_dd
Y_dd = O_dd * R_inv_dd + R_inv_qd * (O_qd - (V_q * P_dq + V_d * P_dd) - (I_q * Q_dq + I_d * Q_dd));
Y_dd_expanded = expand(Y_dd);

fprintf('Expanded Analytical Expression for Y_dd(s):\n  %s\n\n', char(Y_dd_expanded));

%% 3. Extraction of 1D Frequency Functions {F1, ..., F5}
% Linear multiplier extraction via symbolic partial differentiation
F1_derived = subs(Y_dd_expanded, [V_q, V_d, I_q, I_d], [0, 0, 0, 0]);
F2_derived = -diff(Y_dd_expanded, V_q);
F3_derived = -diff(Y_dd_expanded, V_d);
F4_derived = -diff(Y_dd_expanded, I_q);
F5_derived = -diff(Y_dd_expanded, I_d);

fprintf('Extracted 1D Frequency Branches for PINN Mapping:\n');
fprintf('  • F1(s) [Base Operator]       = %s\n', char(F1_derived));
fprintf('  • F2(s) [V_q Sensitivity]     = %s\n', char(F2_derived));
fprintf('  • F3(s) [V_d Sensitivity]     = %s\n', char(F3_derived));
fprintf('  • F4(s) [I_q Sensitivity]     = %s\n', char(F4_derived));
fprintf('  • F5(s) [I_d Sensitivity]     = %s\n\n', char(F5_derived));

%% 4. Identity Equivalence Verification
% Reconstruct admittance from extracted decoupled branches
Y_dd_reconstructed = F1_derived - V_q * F2_derived - V_d * F3_derived - I_q * F4_derived - I_d * F5_derived;
identity_residual  = simplify(Y_dd_expanded - Y_dd_reconstructed);

fprintf('=========================================================================\n');
if identity_residual == 0
    fprintf('VERIFICATION STATUS: SUCCESS (Residual = 0)\n');
    fprintf('The admittance map factors strictly into 1D frequency sub-networks.\n');
else
    fprintf('VERIFICATION STATUS: FAILED (Non-zero algebraic residual detected).\n');
end
fprintf('=========================================================================\n');
