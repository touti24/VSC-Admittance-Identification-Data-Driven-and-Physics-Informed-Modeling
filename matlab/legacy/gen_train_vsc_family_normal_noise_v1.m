%% generate_training_data_vsc_family_normal_noise.m
% =========================================================================
% CORRECTED per review meeting: "Don't try to train your neural network
% under extreme noise... test under normal conditions where the noise
% maximum can reach 3%, 4%, 5%, that's it." Aya's SNR_dB_list previously
% went down to 5dB (~56% relative noise power) -- far past "normal."
% 3-5% relative noise amplitude maps to roughly:
%     SNR_dB = 20*log10(1/relative_noise)
%     3% -> ~30.5 dB,  4% -> ~28 dB,  5% -> ~26 dB
% So the training sweep is corrected to [30 28 26] dB only. Everything
% else (voltage-step global-frame injection, PLL damping fix, physical
% parameter sweep, non-uniform frequency grid, clean copy) is unchanged
% from the validated version.
% =========================================================================

clear; clc;

%% ---------------- CONFIG ----------------
n_systems     = 3000;
Fs            = 2500;
T_window      = 0.98;
step_size_V   = 1.0;
SNR_dB_list   = [30 28 26];   % CORRECTED: normal conditions only (~3-5% noise)
out_dir       = "training_data_out_normal_noise";
rng(7);

Lf_nom = 1e-3;  Rf_nom = 0.3e-3;
Kp_i_nom = 10.5; Ki_i_nom = 5900;
Kp_pll_nom = 1.2; Ki_pll_nom = 257;
range_frac = 0.5;
Id_range   = [0, 100];
Iq_range   = [-50, 50];
Vd_ss = 311;

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
t = (0:1/Fs:T_window-1/Fs)';
Nt = length(t);
w0 = 2*pi*50;
manifest = {};

% non-uniform frequency grid, dense at the resonance band (unchanged)
f_coarse_lo = linspace(1, 35, 60);
f_dense     = linspace(35, 75, 2000);
f_coarse_hi = linspace(75, 100, 60);
freqs_Hz    = unique([f_coarse_lo, f_dense, f_coarse_hi]);

sys_id = 0;
attempts = 0;
max_attempts = n_systems * 3;

while sys_id < n_systems && attempts < max_attempts
    attempts = attempts + 1;

    Lf     = Lf_nom     * (1 + range_frac*(2*rand()-1));
    Rf     = Rf_nom     * (1 + range_frac*(2*rand()-1));
    Kp_i   = Kp_i_nom   * (1 + range_frac*(2*rand()-1));
    Ki_i   = Ki_i_nom   * (1 + range_frac*(2*rand()-1));
    Kp_pll = Kp_pll_nom * (1 + range_frac*(2*rand()-1));
    Ki_pll = Ki_pll_nom * (1 + range_frac*(2*rand()-1));
    Id_eq  = Id_range(1) + diff(Id_range)*rand();
    Iq_eq  = Iq_range(1) + diff(Iq_range)*rand();

    n = 6;
    A = zeros(n); B = zeros(n,2); C = zeros(2,n);
    A(1,1) = -Rf/Lf;   A(1,2) =  w0;
    A(2,1) = -w0;      A(2,2) = -Rf/Lf;
    A(1,3) =  Kp_i/Lf; A(2,4) =  Kp_i/Lf;
    A(3,1) = -1;       A(4,2) = -1;
    A(5,5) = -Kp_pll*Vd_ss;
    A(5,6) = 1;
    A(6,5) = -Ki_pll*Vd_ss;
    A(1,5) =  Iq_eq*w0;
    A(2,5) = -Id_eq*w0;
    B(1,1) = -1/Lf;
    B(2,2) = -1/Lf;
    B(5,2) =  1;
    C(1,1) = 1;
    C(2,2) = 1;
    D = zeros(2,2);

    if any(real(eig(A)) >= 0)
        continue;
    end
    assert(isreal(A) && isreal(B) && isreal(C), "State-space matrices must be real");
    sysL = ss(A,B,C,D);

    sys_id = sys_id + 1;

    u_ev1 = zeros(Nt,2); u_ev1(t>=0.1,1) = step_size_V;
    u_ev2 = zeros(Nt,2); u_ev2(t>=0.1,2) = step_size_V;
    y_ev1_clean = lsim(sysL, u_ev1, t);
    y_ev2_clean = lsim(sysL, u_ev2, t);

    Y_true = zeros(2,2,length(freqs_Hz));
    for k = 1:length(freqs_Hz)
        s = 1i*2*pi*freqs_Hz(k);
        Y_true(:,:,k) = C*((s*eye(n)-A)\B) + D;
    end

    % clean copy first (unchanged from validated version)
    T_clean = table(t, y_ev1_clean(:,1), y_ev1_clean(:,2), y_ev2_clean(:,1), y_ev2_clean(:,2), ...
        'VariableNames', {'time','id_ev1','iq_ev1','id_ev2','iq_ev2'});
    writetable(T_clean, fullfile(out_dir, sprintf('sys_%05d_SNR_Inf.csv', sys_id)));

    for snr_idx = 1:length(SNR_dB_list)
        SNR_dB = SNR_dB_list(snr_idx);
        y_ev1 = add_sensor_noise(y_ev1_clean, t, SNR_dB);
        y_ev2 = add_sensor_noise(y_ev2_clean, t, SNR_dB);

        T = table(t, y_ev1(:,1), y_ev1(:,2), y_ev2(:,1), y_ev2(:,2), ...
            'VariableNames', {'time','id_ev1','iq_ev1','id_ev2','iq_ev2'});
        fname = sprintf('sys_%05d_SNR_%d.csv', sys_id, SNR_dB);
        writetable(T, fullfile(out_dir, fname));
    end

    ymat_fname = sprintf('sys_%05d_Ytrue.mat', sys_id);
    save(fullfile(out_dir, ymat_fname), 'Y_true', 'freqs_Hz', ...
         'Lf','Rf','Kp_i','Ki_i','Kp_pll','Ki_pll','Id_eq','Iq_eq');

    manifest(end+1,:) = {sys_id, Id_eq, Iq_eq, ymat_fname}; %#ok<SAGROW>

    if mod(sys_id, 200) == 0
        fprintf("Generated %d / %d systems\n", sys_id, n_systems);
    end
end

if isempty(manifest)
    error("manifest is EMPTY -- every parameter draw was unstable. Re-check A matrix construction.");
end

writetable(cell2table(manifest, 'VariableNames', {'sys_id','Id_eq','Iq_eq','ytrue_file'}), ...
    fullfile(out_dir, "manifest.csv"));

fprintf("\nDone. %d systems written to '%s/' with NORMAL noise only (SNR = %s dB).\n", ...
    sys_id, out_dir, mat2str(SNR_dB_list));

function y_noisy = add_sensor_noise(y_clean, t, SNR_dB)
    y_noisy = y_clean;
    step_idx = t >= 0.1;
    for ch = 1:size(y_clean,2)
        sig = y_clean(step_idx, ch);
        sig_power = mean((sig - mean(sig)).^2);
        if sig_power == 0
            sig_power = eps;
        end
        noise_power = sig_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power) * randn(size(y_clean,1),1);
        y_noisy(:,ch) = y_clean(:,ch) + noise;
    end
end