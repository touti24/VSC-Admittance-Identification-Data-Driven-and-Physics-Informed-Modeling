%% params_VSC_I.m
% =========================================================================
% VSC I Parameter Initialization Script
%
% Description:
%   Initializes workspace variables for the baseline grid-connected converter
%   (VSC I) using hardware and controller specifications from Table I of
%   Zhang et al., IEEE TIE (2023).
%
% Usage:
%   Run prior to launching Simulink models (VSC-1.slx, VSC_Iavrg.slx)
%   or automated impedance measurement scripts.
% =========================================================================

clc;

%% 1. Grid Specifications (50 Hz Three-Phase System)
f0                        = 50;                     % Nominal grid frequency [Hz]
omega0                    = 2 * pi * f0;            % Angular frequency [rad/s]
V_grid_rms                = 380;                    % Line-to-line RMS voltage [V]
V_grid_phase_to_phase_rms = 380;                    % Reference line-to-line RMS [V]
V_grid_peak               = V_grid_rms * sqrt(2/3); % Peak phase-to-neutral voltage [V]

%% 2. Power Stage & Hardware Parameters
Vdc                       = 700;                    % DC-link voltage [V]
fsw                       = 20000;                  % Inverter switching frequency [Hz]
Ts                        = 1 / (10 * fsw);         % Simulation integration timestep [s] (5e-6 s)

%% 3. AC Filter Specification (L-Filter)
Lf                        = 1.0e-3;                 % Filter inductance [H]
Rf                        = 0.3e-3;                 % Parasitic resistance [Ohm]

%% 4. Controller Parameters (Synchronous dq-Frame)
% Inner-loop PI current controller
Kp_i                      = 10.5;                   % Proportional gain
Ki_i                      = 5900;                   % Integral gain [rad/s]

% Synchronous Reference Frame Phase-Locked Loop (SRF-PLL)
Kp_pll                    = 1.2;                    % Proportional gain [rad/(V*s)]
Ki_pll                    = 257;                    % Integral gain [rad/(V*s^2)]

%% 5. Default Steady-State Operating Point
Id_ref                    = 40;                     % Active current reference [A]
Iq_ref                    = 0;                      % Reactive current reference [A] (unity PF)

%% 6. Small-Signal Perturbation Defaults
% Fallback definitions when not executing an automated sweep
if ~exist('f_p', 'var')
    fp                    = 10;                     % Perturbation frequency [Hz]
    perturb_amp           = 0;                      % Perturbation amplitude [A]
end

fprintf('[params_VSC_I] Workspace successfully initialized with VSC I parameters.\n');
