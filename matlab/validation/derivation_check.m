% ==============================================================================
% phase4_pinn/derivation_check.m
% Verification of Admittance Factorization Architecture (Eqs. 13-18)
% ==============================================================================
clear;
clc;

disp('--- Running MATLAB Symbolic Physics Factorization Verification ---');
disp(' ');

% 1. Define independent physical operating point scalars
syms f positive
syms V_d V_q I_d I_q real

% 2. Define the individual frequency-dependent variables from Eq. 15
syms O_dd O_dq O_qd O_qq real
syms P_dd P_dq P_qd P_qq real
syms Q_dd Q_dq Q_qd Q_qq real
syms R_inv_dd R_inv_dq R_inv_qd R_inv_qq real

% 3. Construct the analytical algebraic combiner equation for Y_dd based on Eq. 16
% Y_dd = O_dd*R_inv_dd + R_inv_qd * (O_qd - (V_q*P_dq + V_d*P_dd) - (I_q*Q_dq + I_d*Q_dd))
Y_dd = O_dd*R_inv_dd + R_inv_qd * (O_qd - (V_q*P_dq + V_d*P_dd) - (I_q*Q_dq + I_d*Q_dd));

% Expand the expression to isolate multi-variable products
Y_dd_expanded = expand(Y_dd);

disp('Expanded Equation for Y_dd(s):');
disp(Y_dd_expanded);
disp('----------------------------------------------------------------------');

% 4. Differentiate to verify the linear multipliers (F1 to F5 coefficients)
% Y_dd = F1 - V_q*F2 - V_d*F3 - I_q*F4 - I_d*F5
F1_derived = subs(Y_dd_expanded, [V_q, V_d, I_q, I_d], [0, 0, 0, 0]);
F2_derived = -diff(Y_dd_expanded, V_q);
F3_derived = -diff(Y_dd_expanded, V_d);
F4_derived = -diff(Y_dd_expanded, I_q);
F5_derived = -diff(Y_dd_expanded, I_d);

disp('Extracted Frequency Functions for PINN Mapping:');
fprintf('  • F1(s) [Base Term]     = %s\n', char(F1_derived));
fprintf('  • F2(s) [V_q Coefficient] = %s\n', char(F2_derived));
fprintf('  • F3(s) [V_d Coefficient] = %s\n', char(F3_derived));
fprintf('  • F4(s) [I_q Coefficient] = %s\n', char(F4_derived));
fprintf('  • F5(s) [I_d Coefficient] = %s\n', char(F5_derived));
disp('----------------------------------------------------------------------');

% 5. Reconstruct and perform mathematical identity check
Y_dd_reconstructed = F1_derived - V_q*F2_derived - V_d*F3_derived - I_q*F4_derived - I_d*F5_derived;
identity_check = simplify(Y_dd_expanded - Y_dd_reconstructed);

if identity_check == 0
    disp('✅ VERIFICATION SUCCESSFUL: Equation 17 holds perfectly.');
    disp('The system reduces cleanly to 1D frequency functions.');
else
    disp('❌ Verification mismatch detected.');
end