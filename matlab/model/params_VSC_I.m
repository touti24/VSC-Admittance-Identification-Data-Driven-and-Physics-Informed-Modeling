%% params_VSC_I.m
% Sets up the workspace variables for VSC I (Base simulation model)
clc; 

%% Grid Configuration (50 Hz System)
f0 = 50;                % Fundamental frequency (Hz)
omega0 = 2*pi*f0;       % Fundamental angular frequency (rad/s)
V_grid_rms = 380;       % Standard Line-to-Line RMS voltage (V)
V_grid_phase_to_phase_rms = 380; 
V_grid_peak = V_grid_rms * sqrt(2/3); 

%% Hardware & Power Stage Constants
Vdc = 700;              % DC link voltage (V)
fsw = 20000;            % Switching frequency (20 kHz)
Ts = 1 / (10 * fsw);    % Ts = 5e-6 s

%% L-Filter Specification (Table I)
Lf = 1e-3;              % 1 mH Inverter inductor
Rf = 0.3e-3;            % 0.3 mOhm parasitic resistance

%% Inner Current Loop Controllers (Table I)
Kp_i = 10.5;            % Proportional gain
Ki_i = 5900;            % Integral gain

%% Phase-Locked Loop (SRF-PLL) Controllers (Table I)
Kp_pll = 1.2;           % PLL Proportional gain
Ki_pll = 257;           % PLL Integral gain

%% Steady-State Operating Point (Defaults)
Id_ref = 40;            % Start tracking at 40 A active current
Iq_ref = 0;             % Unity power factor

%% Dynamic Injection Channels (Only set if not running a sweep)
if ~exist('f_p', 'var')
    fp = 10;                
    perturb_amp = 0;        
end
fprintf('Workspace initialized with VSC I parameters.\n');