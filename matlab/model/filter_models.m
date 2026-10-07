function [Yp, Yo] = filter_models(filter_type, params)
% =========================================================================
% filter_models.m
%
% Returns the filter admittance transfer functions Yp(s) and Yo(s)
% as MATLAB tf objects, for use in the analytical admittance model.
%
% The paper uses two filter topologies:
%   'L'   — single inductor (VSC I and VSC II)
%   'LCL' — LCL filter with capacitor-current active damping (Phase 6)
%
% ROLE IN EQ. (12):
%   Yp(s) = admittance seen from the converter output (inverter side)
%   Yo(s) = admittance seen from the PCC (grid side)
%   For L-filter: Yp = Yo = 1/(Lf*s + Rf)
%   For LCL:      Yp = grid-side admittance, Yo = converter-side admittance
%                 (more complex — see LCL section below)
%
% INPUTS:
%   filter_type — 'L' or 'LCL'  (string)
%   params      — struct:
%     For 'L':
%       .Lf   [H]    inverter inductor
%       .Rf   [Ohm]  parasitic resistance
%     For 'LCL' (additional fields):
%       .L2   [H]    grid-side inductor
%       .Cf   [F]    filter capacitor
%       .Kd         active damping gain (capacitor current feedback)
%
% OUTPUTS:
%   Yp  — MATLAB tf object  (converter-side admittance)
%   Yo  — MATLAB tf object  (grid-side / PCC admittance)
%
% USAGE:
%   params = struct('Lf',1e-3,'Rf',3e-4);
%   [Yp, Yo] = filter_models('L', params);
%   bode(Yp);
% =========================================================================

    switch lower(filter_type)

        % ── L-FILTER (VSC I and VSC II) ───────────────────────────────────
        case 'l'
            Lf = params.Lf;
            Rf = params.Rf;

            % Yp(s) = 1 / (Lf*s + Rf)
            Yp = tf(1, [Lf, Rf]);

            % For L-filter: Yo = Yp (single branch, no capacitor)
            Yo = Yp;

            fprintf('filter_models: L-filter | Lf=%.2f mH | Rf=%.3f mOhm\n', ...
                    Lf*1e3, Rf*1e3);

        % ── LCL-FILTER WITH ACTIVE DAMPING (Phase 6) ──────────────────────
        % Structure: Converter → L1(=Lf,Rf) → Capacitor Cf → L2 → Grid
        % Active damping: capacitor current ×Kd fed back to subtract from
        % modulation voltage (reduces the resonance peak without dissipation)
        case 'lcl'
            Lf = params.Lf;     % converter-side inductor [H]
            Rf = params.Rf;     % converter-side parasitic [Ohm]
            L2 = params.L2;     % grid-side inductor [H]   (= 3 mH, paper §IV-B)
            Cf = params.Cf;     % filter capacitor [F]     (= 5 µF, paper §IV-B)
            Kd = params.Kd;     % active damping gain      (= 10, assumption)

            % Without active damping, LCL admittance (converter-side):
            %   Yp_nodamp = Cf*s / (Lf*Cf*s^2 + Rf*Cf*s + 1) ... complex
            %
            % With capacitor-current active damping, the effective admittance
            % is modified. For the impedance measurement purposes we use
            % the output admittance (grid-side):
            %
            %   Yo(s) = 1 / (L2*s + Lf*s + Rf + Kd)
            %   (simplified: active damping acts like an added series resistance Kd)
            %
            % NOTE: This is an approximation. For full derivation see
            % Fig. 16 of the paper and Ref [16] (Pan et al. 2015).

            % Converter-side admittance (simplified with active damping)
            % Yp(s) = 1/(Lf*s + Rf + Kd)
            Yp = tf(1, [Lf, Rf + Kd]);

            % Grid-side admittance
            % Yo(s) = 1/(L2*s)  ... grid-side inductor only (no damping on grid side)
            Yo = tf(1, [L2, 0]);

            fprintf('filter_models: LCL-filter | Lf=%.0f mH | L2=%.0f mH | Cf=%.0f uF | Kd=%.0f\n', ...
                    Lf*1e3, L2*1e3, Cf*1e6, Kd);
            fprintf('  NOTE: Active damping modelled as series resistance Kd=%.0f Ohm (approximation)\n', Kd);

        otherwise
            error('filter_models: unknown filter type ''%s''. Use ''L'' or ''LCL''.', filter_type);
    end
end