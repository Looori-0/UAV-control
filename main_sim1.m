function main_sim1()
clc; clear; close all;

%% ================== Sim-1: no disturbance + large initial error ==================
sim.Tend = 8;
sim.dt   = 0.001;

% inertia
sim.J = diag([0.054, 0.053, 0.089]);

% desired attitude (constant)
sim.qd = [1;0;0;0];

% large initial error via axis-angle (avoid Euler singularity)
axis0  = [1; -0.3; 0.5];
axis0  = axis0 / norm(axis0);
theta0 = deg2rad(100);               % large error
sim.q0  = [cos(theta0/2); axis0*sin(theta0/2)];

% initial angular rate (non-zero)
sim.w0 = [0.2; -0.15; 0.10];

% torque saturation (engineering constraint)
sim.tau_limit = [8; 8; 8];           % N·m

% disturbance function (Sim-1: d=0)
sim.dist_fn = @(t,x) [0;0;0];

%% ================== PID parameters ==================
sim.pid.Kp = diag([6, 6, 4]);
sim.pid.Kd = diag([2.0, 2.0, 1.5]);
sim.pid.Ki = diag([1.0, 1.0, 0.8]);
sim.pid.int_limit = [0.5; 0.5; 0.5];

%% ================== SMC parameters ==================
sim.smc.c   = 4.0;
sim.smc.K   = diag([3.0, 3.0, 2.5]);
sim.smc.eta = diag([0.6, 0.6, 0.5]);
sim.smc.phi = 0.05;                  % boundary layer

%% ================== TSMC parameters (baseline) ==================
% terminal surface: s = e_dot + lambda * sig(e)^r, 0<r<1
% reaching: s_dot = -K * sat(s/phi)
sim.tsmc.lambda = 3.0;
sim.tsmc.r      = 0.6;
sim.tsmc.K      = diag([4.0, 4.0, 3.5]);
sim.tsmc.phi    = 0.05;
% numerical fallback when A is ill-conditioned
sim.tsmc.mu        = 1e-4;
sim.tsmc.rcond_thr = 1e-8;

%% ================== FFTSMC parameters ==================
% s = e_dot + lambda1*sig(e)^(p/q) + lambda2*sig(e)^(q/p)
sim.ffts.p = 1; 
sim.ffts.q = 2;                       % r1=0.5, r2=2

sim.ffts.lambda1 = 3.0;
sim.ffts.lambda2 = 1.2;

% reaching law: -K1*sig(s)^alpha - K2*sig(s)^beta
sim.ffts.alpha = 0.6;                 % 0<alpha<1
sim.ffts.beta  = 1.2;                 % beta>1
sim.ffts.K1 = diag([6.0, 6.0, 5.0]);
sim.ffts.K2 = diag([1.5, 1.5, 1.2]);

% DLS fallback when A is ill-conditioned
sim.ffts.mu = 1e-4;
sim.ffts.rcond_thr = 1e-8;

%% ================== FTDO parameters (Eq. 3.1) ==================
% tilde = w - what
% what_dot = J^{-1}(tau - w×(Jw) + dhat) + k1*sqrt(|tilde|)*sgn(tilde)
% dhat_dot = k2*sgn(tilde)
% (use sat(tilde/phi) instead of pure sign for numerical stability)
sim.ftdo.phi = 0.03;
sim.ftdo.k1  = diag([6, 6, 6]);
sim.ftdo.k2  = diag([10,10,10]);

%% ================== Run ==================
res_pid  = attitude_sim('PID', sim);
res_smc  = attitude_sim('SMC', sim);
res_tsmc = attitude_sim('TSMC', sim);
res_ffts = attitude_sim('FFTSMC_FTDO', sim);

t = res_pid.t;

%% ================== Figure 1: attitude error angle theta_e(t) ==================
figure('Name','Attitude Error Angle','Color','w');
plot(t, rad2deg(res_ffts.theta_e),'LineWidth',1.3); hold on;
plot(t, rad2deg(res_tsmc.theta_e),'LineWidth',1.3);
plot(t, rad2deg(res_pid.theta_e),'LineWidth',1.3);
plot(t, rad2deg(res_smc.theta_e),'LineWidth',1.3);
grid on; xlabel('Time (s)'); ylabel('\theta_e (deg)');
title('Attitude Error Angle \theta_e(t)');
legend('FFTSMC+FTDO','TSMC','PID','SMC','Location','best');

%% ================== Figure 2: error quaternion vector part e=[e1,e2,e3] ==================
figure('Name','Quaternion Error Vector Part','Color','w');
for i=1:3
    subplot(3,1,i);
    plot(t, res_ffts.ev(:,i),'LineWidth',1.2); hold on;
    plot(t, res_tsmc.ev(:,i),'LineWidth',1.2);
    plot(t, res_pid.ev(:,i),'LineWidth',1.2);
    plot(t, res_smc.ev(:,i),'LineWidth',1.2);
    grid on; ylabel(sprintf('e_%d',i));
    if i==1, title('Error Quaternion Vector Part e(t)'); end
    if i==3, xlabel('Time (s)'); end
end
legend('FFTSMC+FTDO','TSMC','PID','SMC','Location','best');

%% ================== Figure 3: angular rate omega ==================
figure('Name','Angular Rate','Color','w');
labels = {'\\omega_x','\\omega_y','\\omega_z'};
for i=1:3
    subplot(3,1,i);
    plot(t, res_ffts.w(:,i),'LineWidth',1.2); hold on;
    plot(t, res_tsmc.w(:,i),'LineWidth',1.2);
    plot(t, res_pid.w(:,i),'LineWidth',1.2);
    plot(t, res_smc.w(:,i),'LineWidth',1.2);
    grid on; ylabel([labels{i} ' (rad/s)']);
    if i==1, title('Angular Rate \\omega(t)'); end
    if i==3, xlabel('Time (s)'); end
end
legend('FFTSMC+FTDO','TSMC','PID','SMC','Location','best');

%% ================== Figure 4: control torque tau ==================
figure('Name','Control Torque','Color','w');
labels = {'\\tau_x','\\tau_y','\\tau_z'};
for i=1:3
    subplot(3,1,i);
    plot(t, res_ffts.tau(:,i),'LineWidth',1.2); hold on;
    plot(t, res_tsmc.tau(:,i),'LineWidth',1.2);
    plot(t, res_pid.tau(:,i),'LineWidth',1.2);
    plot(t, res_smc.tau(:,i),'LineWidth',1.2);
    grid on; ylabel([labels{i} ' (N\\cdotm)']);
    if i==1, title('Control Torque \\tau(t)'); end
    if i==3, xlabel('Time (s)'); end
end
legend('FFTSMC+FTDO','TSMC','PID','SMC','Location','best');

%% ================== quick summary ==================
fprintf('Final theta_e (deg): FFTS=%.4f, TSMC=%.4f, PID=%.4f, SMC=%.4f\n', ...
    rad2deg(res_ffts.theta_e(end)), rad2deg(res_tsmc.theta_e(end)), rad2deg(res_pid.theta_e(end)), rad2deg(res_smc.theta_e(end)));

end
