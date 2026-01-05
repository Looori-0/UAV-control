function out = attitude_sim(ctrlName, sim)
% attitude_sim - Quaternion attitude simulation with PID / SMC / TSMC / FFTSMC+FTDO
%
% ctrlName: 'PID' | 'SMC' | 'TSMC' | 'FFTSMC_FTDO'
% sim (struct) required fields:
%   Tend, dt, J(3x3), qd(4x1), q0(4x1), w0(3x1)
%   pid:  Kp,Kd,Ki,int_limit
%   smc:  c,K,eta,phi
%   tsmc: lambda,r,K,phi,mu,(optional) rcond_thr
%   ffts: p,q,lambda1,lambda2,alpha,beta,K1,K2,mu,(optional) rcond_thr
%   ftdo: k1,k2,phi   (Eq.(3.1) super-twisting form; phi used for sat instead of pure sign)
% optional:
%   tau_limit (scalar or 3x1), dist_fn (function handle d(t,x)->3x1)
%
% Output:
%   out.t, out.q, out.w, out.tau, out.ev, out.theta_e
%   out.dhat (only FFTSMC_FTDO)

    % ---------------- defaults ----------------
    if ~isfield(sim,'tau_limit'), sim.tau_limit = inf; end
    if ~isfield(sim,'dist_fn') || isempty(sim.dist_fn)
        sim.dist_fn = @(t,x) [0;0;0]; % default: no disturbance
    end
    if ~isfield(sim,'ffts') || ~isfield(sim.ffts,'rcond_thr')
        sim.ffts.rcond_thr = 1e-8;
    end
    if ~isfield(sim,'tsmc') || ~isfield(sim.tsmc,'rcond_thr')
        sim.tsmc.rcond_thr = 1e-8;
    end

    Tend = sim.Tend;
    dt   = sim.dt;

    qd = sim.qd(:);
    q0 = quat_normalize(sim.q0(:));
    w0 = sim.w0(:);

    t_grid = (0:dt:Tend)';

    % ---------------- initial state ----------------
    switch upper(ctrlName)
        case 'PID'
            x0 = [q0; w0; zeros(3,1)];          % [q(4); w(3); int_e(3)]
        case {'SMC','TSMC'}
            x0 = [q0; w0];                       % [q(4); w(3)]
        case 'FFTSMC_FTDO'
            what0 = w0;                          % hat(omega)
            dhat0 = zeros(3,1);                  % hat(d)
            x0 = [q0; w0; what0; dhat0];         % [q; w; what; dhat]
        otherwise
            error('Unknown ctrlName: %s', ctrlName);
    end

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

    % ---------------- integrate ----------------
    [t_ode, X_ode] = ode45(@(t,x) dyn_all(t, x, sim, ctrlName), [0 Tend], x0, opts);

    % ---------------- interpolate to uniform grid ----------------
    Xg = interp1(t_ode, X_ode, t_grid, 'pchip');

    % ---------------- parse outputs ----------------
    q = normalize_rows_quat(Xg(:,1:4));
    w = Xg(:,5:7);

    [~, ev, theta_e] = quat_error_series(q, qd);

    % ---------------- compute tau (use ODE nodes then interpolate) ----------------
    N_ode = size(X_ode,1);
    tau_ode = zeros(N_ode,3);
    for i = 1:N_ode
        tau_i = control_tau(ctrlName, X_ode(i,:).', sim);
        tau_ode(i,:) = tau_i(:).';
    end
    % interpolate torque to uniform grid (linear avoids artifacts for non-smooth control)
    tau = interp1(t_ode, tau_ode, t_grid, 'linear');


    out.t = t_grid;
    out.q = q;
    out.w = w;
    out.tau = tau;
    out.ev = ev;
    out.theta_e = theta_e;

    if strcmpi(ctrlName,'FFTSMC_FTDO')
        out.dhat = interp1(t_ode, X_ode(:,11:13), t_grid, 'linear');
    end
end

% ======================================================================
%                           Dynamics
% ======================================================================
function dx = dyn_all(t, x, sim, ctrlName)
    J  = sim.J;
    qd = sim.qd(:);

    q = quat_normalize(x(1:4));
    w = x(5:7);

    d = sim.dist_fn(t, x);   % 3x1 disturbance torque

    tau = control_tau(ctrlName, x, sim);

    % rigid body: J*w_dot = tau - w×(Jw) + d
    w_dot = J \ (tau - cross(w, J*w) + d);

    % quaternion kinematics
    q_dot = quat_kinematics(q, w);

    switch upper(ctrlName)
        case 'PID'
            qe = quat_mul(quat_conj(qd), q);
            qe = shortest_quat(qe);
            e  = qe(2:4);

            int_lim = sim.pid.int_limit(:);
            int_e_dot = clamp_vec(e, int_lim);

            dx = [q_dot; w_dot; int_e_dot];

        case {'SMC','TSMC'}
            dx = [q_dot; w_dot];

        case 'FFTSMC_FTDO'
            % ---------------- FTDO (Eq. 3.1) ----------------
            % tilde = w - what
            what = x(8:10);
            dhat = x(11:13);

            tilde = w - what;

            % sgn(tilde) implemented with sat boundary layer for numerical stability
            phi_o = sim.ftdo.phi;
            sgn_tilde = sat_vec(tilde, phi_o);

            k1 = gain3(sim.ftdo.k1);
            k2 = gain3(sim.ftdo.k2);

            what_dot = (J \ (tau - cross(w, J*w) + dhat)) ...
                     + k1 * (sqrt(abs(tilde) + 1e-12) .* sgn_tilde);

            dhat_dot = k2 * sgn_tilde;

            dx = [q_dot; w_dot; what_dot; dhat_dot];

        otherwise
            error('Unknown ctrlName in dyn_all');
    end
end

% ======================================================================
%                           Control
% ======================================================================
function tau = control_tau(ctrlName, x, sim)
    J  = sim.J;
    qd = sim.qd(:);

    q = quat_normalize(x(1:4));
    w = x(5:7);

    % error quaternion qe = qd^{-1} ⊗ q
    qe = quat_mul(quat_conj(qd), q);
    qe = shortest_quat(qe);

    qe0 = qe(1);
    e   = qe(2:4);

    % e_dot = A(qe) * w
    A = 0.5*(qe0*eye(3) + skew3(e));
    e_dot = A*w;

    switch upper(ctrlName)
        case 'PID'
            int_e = x(8:10);
            Kp = gain3(sim.pid.Kp);
            Kd = gain3(sim.pid.Kd);
            Ki = gain3(sim.pid.Ki);
            tau = -Kp*e - Kd*w - Ki*int_e;

        case 'SMC'
            c   = sim.smc.c;
            K   = gain3(sim.smc.K);
            eta = gain3(sim.smc.eta);
            phi = sim.smc.phi;

            s = w + c*e;
            sat_s = sat_vec(s, phi);

            % enforce: w_dot = -c*e_dot - K*s - eta*sat(s)
            w_dot_cmd = -c*e_dot - K*s - eta*sat_s;

            tau = cross(w, J*w) + J*w_dot_cmd;

        case 'TSMC'
            % ---------------- Terminal Sliding Mode (baseline) ----------------
            % terminal sliding surface: s = e_dot + lambda * sig(e)^r, 0<r<1
            lam = sim.tsmc.lambda;
            r   = sim.tsmc.r;          % 0<r<1
            K   = gain3(sim.tsmc.K);
            phi = sim.tsmc.phi;

            s = e_dot + lam * sig_pow(e, r);

            % G_t(e) = diag(lam*r*|e|^(r-1))
            eps0 = 1e-8;
            Gdiag = lam*r*(abs(e)+eps0).^(r-1);

            % Adot = 0.5*(qe0_dot*I + skew(e_dot)), qe0_dot = -0.5*e^T*w
            qe0_dot = -0.5*(e.'*w);
            Adot = 0.5*(qe0_dot*eye(3) + skew3(e_dot));

            % reaching: s_dot = -K*sat(s/phi)
            reach = K * sat_vec(s, phi);

            % v = -Adot*w - G(e)*e_dot - reach
            v = -Adot*w - (Gdiag .* e_dot) - reach;

            % solve A*w_dot_cmd = v  (standard A^{-1}v, with smooth DLS blending)
            mu   = sim.tsmc.mu;
            thr1 = sim.tsmc.rcond_thr;
            thr2 = 10*thr1;

            w_dot_cmd = A \ v;
            w2 = (A.'*A + mu*eye(3)) \ (A.'*v); % DLS

            if any(~isfinite(w_dot_cmd))
                w_dot_cmd = w2;
            end

            cnd = rcond(A);
            g = (cnd - thr1) / max(thr2 - thr1, 1e-12);
            g = min(max(g, 0), 1);

            w_dot_cmd = g*w_dot_cmd + (1-g)*w2;

            tau = cross(w, J*w) + J*w_dot_cmd;

        case 'FFTSMC_FTDO'
      
           % ---------------- FFTSMC ----------------
            p  = sim.ffts.p;
            qn = sim.ffts.q;
            r1 = p/qn;     % 0.5
            r2 = qn/p;     % 2.0

            lam1 = sim.ffts.lambda1;
            lam2 = sim.ffts.lambda2;

            % sliding surface (保持不变)
            s = e_dot + lam1*sig_pow(e, r1) + lam2*sig_pow(e, r2);

            eps_reg = 0.002;  % 正则化因子，数值越大曲线越平滑，建议 0.002 ~ 0.01
            
            % 第一项 (奇异项): r1 * (|e| + eps_reg)^(r1-1)
            term1 = r1 * (abs(e) + eps_reg).^(r1-1);
            
            % 第二项 (非奇异项): r2 * |e|^(r2-1) (r2=2, r2-1=1，本身就是线性的，无需正则)
            term2 = r2 * (abs(e)).^(r2-1);
            
            Gdiag = lam1*term1 + lam2*term2;

            % Adot = 0.5*(qe0_dot*I + skew(e_dot))
            qe0_dot = -0.5*(e.'*w);
            Adot = 0.5*(qe0_dot*eye(3) + skew3(e_dot));

            % reaching law: -K1 sig(s)^alpha - K2 sig(s)^beta
            alpha = sim.ffts.alpha;
            beta  = sim.ffts.beta;
            K1 = gain3(sim.ffts.K1);
            K2 = gain3(sim.ffts.K2);

            reach = K1*sig_pow(s, alpha) + K2*sig_pow(s, beta);

            % v = -Adot*w - G(e)*e_dot - reach
            v = -Adot*w - (Gdiag .* e_dot) - reach;

            mu   = sim.ffts.mu;
            thr1 = sim.ffts.rcond_thr;
            thr2 = 10*thr1;

            w_dot_cmd = A \ v;

            % use dhat from observer state
            dhat = x(11:13);

            tau = cross(w, J*w) - dhat + J*w_dot_cmd;

        otherwise
            error('Unknown ctrlName in control_tau');
    end

    % torque saturation
    tau = clamp_vec(tau, sim.tau_limit);
end

% ======================================================================
%                         Utilities
% ======================================================================
function K = gain3(Kin)
    % Accept scalar, 3x1, or 3x3
    if isscalar(Kin)
        K = Kin*eye(3);
    else
        Kin = Kin(:,:);
        if all(size(Kin)==[3 1])
            K = diag(Kin(:));
        elseif all(size(Kin)==[1 3])
            K = diag(Kin(:).');
        else
            K = Kin;
        end
    end
end

function qdot = quat_kinematics(q, w)
    q0 = q(1); qv = q(2:4);
    qdot0 = -0.5 * (qv.' * w);
    qdotv =  0.5 * (q0*eye(3) + skew3(qv)) * w;
    qdot  = [qdot0; qdotv];
end

function S = skew3(v)
    S = [  0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2)  v(1)   0  ];
end

function qn = quat_normalize(q)
    qn = q ./ max(norm(q), 1e-12);
end

function Qn = normalize_rows_quat(Q)
    Qn = Q;
    for i = 1:size(Q,1)
        Qn(i,:) = Q(i,:) ./ max(norm(Q(i,:)), 1e-12);
    end
end

function qcon = quat_conj(q)
    qcon = [q(1); -q(2:4)];
end

function qc = quat_mul(a, b)
    a0=a(1); av=a(2:4);
    b0=b(1); bv=b(2:4);
    c0 = a0*b0 - dot(av,bv);
    cv = a0*bv + b0*av + cross(av,bv);
    qc = [c0; cv];
end

function qe = shortest_quat(qe)
    if qe(1) < 0
        qe = -qe;
    end
end

function y = sig_pow(x, r)
    % sig(x)^r = |x|^r * sign(x)
    y = (abs(x)).^r .* sign(x);
end

function y = sat_vec(x, phi)
    % sat(x/phi) in [-1,1] elementwise
    z = x ./ max(phi, 1e-12);
    y = max(-1, min(1, z));
end

function y = clamp_vec(x, lim)
    if isscalar(lim)
        y = min(max(x, -lim), lim);
    else
        lim = lim(:);
        y = min(max(x, -lim), lim);
    end
end

function [QE, EV, TH] = quat_error_series(Q, qd)
    N = size(Q,1);
    QE = zeros(N,4);
    EV = zeros(N,3);
    TH = zeros(N,1);

    qd_conj = quat_conj(qd(:));
    for i=1:N
        q = Q(i,:)';
        qe = quat_mul(qd_conj, q);
        qe = shortest_quat(qe);

        QE(i,:) = qe.';
        EV(i,:) = qe(2:4).';

        q0 = max(-1,min(1,qe(1)));
        qv = qe(2:4);
        TH(i) = 2*atan2(norm(qv), q0);
    end
end
