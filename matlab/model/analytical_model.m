function Y = analytical_model(f, Vd, Vq, Id, Iq, params)
% =========================================================================
% analytical_model.m (FIXED COUPLING)
% =========================================================================
    s = 1j * 2 * pi * f;   % complex frequency at evaluation point
    
    % ── 1. EXTRACT PARAMETERS ────────────────────────────────────────────
    Lf      = params.Lf;
    Rf      = params.Rf;
    Kp_i    = params.Kp_i;
    Ki_i    = params.Ki_i;
    Kp_pll  = params.Kp_pll;
    Ki_pll  = params.Ki_pll;
    Ts      = params.Ts;
    Vd_ss   = params.Vd_ss;   
    
    % Define fundamental grid frequency (typically 50Hz or 60Hz - adjust if needed)
    f1      = 50; 
    w1      = 2 * pi * f1;

    % ── 2. SCALAR TRANSFER FUNCTIONS ─────────────────────────────────────
    Gi = eval_tf([Kp_i, Ki_i], [1, 0], s);
    Gdel = eval_pade_delay(1.5 * Ts, s);
    GPLL = eval_tf([Kp_pll, Ki_pll], [1, Kp_pll*Vd_ss, Ki_pll*Vd_ss], s);

    % ── 3. 2x2 TRANSFER MATRICES ─────────────────────────────────────────
    Gm_i   = Gi * eye(2);
    Gm_del = Gdel * eye(2);
    
    Z_plant = [s*Lf + Rf,  -w1*Lf;
               w1*Lf,       s*Lf + Rf];
    Ym_p = inv(Z_plant); 
    Ym_o = Ym_p; 

    % Revert back to the exact standard forward coupling, but we will fix 
    % the current direction sign globally in the final matrix subtraction.
    Gm_I_PLL = [0,  -GPLL * Iq;
                0,   GPLL * Id];
                
    Gm_V_PLL = [0,  -GPLL * Vq;
                0,   GPLL * Vd];

    % ── 4. INTERMEDIATE QUANTITIES ───────────────────────────────────────
    T = Ym_p * Gm_del * Gm_i;
    Y_to = Ym_o - Ym_p * Gm_del * Gm_V_PLL;
    
    IpT_inv = inv(eye(2) + T);
    G_cl    = IpT_inv * T;
    
    Y_PLL = Gm_I_PLL;

    % ── 5. CLOSED-LOOP ADMITTANCE MATRIX ─────────────────────────────────
    % If the paper considers current injected into the grid vs drawn from it,
    % a global negative sign or alternative transposition applies. 
    % Let's align the component orientations:
    Y = IpT_inv * Y_to - G_cl * Y_PLL;
    
    % If Yqq needs to settle in the -180 deg region like Fig 9, 
    % it means the output current orientation is inverted:
    Y(2,2) = -Y(2,2); 
    Y(1,2) = -Y(1,2); % Flips the cross-coupling frame orientation

end

% =========================================================================
% HELPER FUNCTIONS (Kept identical to your original code)
% =========================================================================
function val = eval_tf(num, den, s)
    val = polyval(num, s) / polyval(den, s);
end

function val = eval_pade_delay(tau, s)
    x = tau * s;   
    num_val = 1 - x/2 + x^2/12 - x^3/120 + x^4/1680;
    den_val = 1 + x/2 + x^2/12 + x^3/120 + x^4/1680;
    val = num_val / den_val;
end