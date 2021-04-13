clc; clear; close all;

g = 9.81;
m = 1;
zc = 1;

%%
% q0 = [-0.5, 2]';
% Ts = 1;
% Td = 0.5;
% dt = 0.01;

q0 = [-0.5, 2]';
Ts = 1;
Td = 0.2;
dt = 0.01;

p = {g, m, zc, q0, Ts, Td};

[t_s, q_s] = ode45(@(t, q) LIP_single(t, q, p), 0:dt:Ts, q0);
[t_d, q_d] = ode45(@(t, q) LIP_double(t, q, p), Ts:dt:(Ts+Td), q_s(end,:)');

figure(1);
subplot(2,1,1); hold on;
plot(t_s, q_s(:,1), 'b', t_d, q_d(:,1), 'r');
xlabel('t'); ylabel('x');
subplot(2,1,2); hold on;
plot(t_s, q_s(:,2), 'b', t_d, q_d(:,2), 'r');
xlabel('t'); ylabel('Vx');

tt = 0:dt:(Ts+Td);
q_traj = nan(length(tt),2);
for i = 1:length(tt)
    q_traj(i,:) = gen_traj(q0, Ts, Td, tt(i));
end
subplot(2,1,1), plot(tt, q_traj(:,1), '--k');
subplot(2,1,2), plot(tt, q_traj(:,2), '--k');
legend('single stance', 'double stance', 'trajectory')

%%
function dq = LIP_single(t, q, p)
    % q = [x, x_dot]'
    % u = [u1, u2]'
    % p = {g, m, zc, q0, Ts, Td};

    [g, m, zc, q0, Ts, Td] = deal(p{:});
    kp = 300; kd = 100;
    
    traj = gen_traj(q0, Ts, Td, t);
    u(1) = [kp, kd] * (traj - q);
    u(2) = 0;
    
    dq(1,1) = q(2);
    dq(2,1) = g/zc*q(1) + 1/m/zc*u(1) + 1/m/zc*u(2);
end

function dq = LIP_double(t, q, p)
    % q = [x, x_dot]'
    % u = [u1, u2]'
    % p = {g, m, zc, q0, Ts, Td};

    [g, m, zc, q0, Ts, Td] = deal(p{:});
    kp = 300; kd = 100;
    
    traj = gen_traj(q0, Ts, Td, t);
    uu = [kp, kd] * (traj - q);
    u(1) = (1-(t-Ts)/Td) * uu;
    u(2) = (t-Ts)/Td * uu;
    
    dq(1,1) = q(2);
    dq(2,1) = g/zc*q(1) + 1/m/zc*u(1) + 1/m/zc*u(2);
end

%% backup
% function dq = LIP_func(t, q, u, p)
%     % q = [x, x_dot, isdouble]'
%     % u = [u1, u2]'
%     % p = [g, m, zc];
%     g = p(1);
%     m = p(2);
%     zc = p(3);
%     
%     dq(1,1) = q(2);
%     dq(2,1) = g/zc*q(1) + 1/m/zc*u(1) + 1/m/zc*u(2);
% end

% function dq = LIP_func(t, q, p)
%     % q = [x, x_dot]'
%     % u = [u1, u2]'
%     % p = {g, m, zc, q0, Ts, Td};
% 
%     [g, m, zc, q0, Ts, Td] = deal(p{:});
%     kp = 1; kd = 0.01;
%     
%     t = mod(t, Ts+Td);
%     traj = gen_traj(q0, Ts, Td, t);
%     if t>=0 && t<Ts
%         u(1) = [kp, kd] * (q - traj);
%         u(2) = 0;
%         
%     elseif t>=Ts && t<=(Ts+Td)
%         uu = [kp, kd] * (q - traj);
%         u(1) = (1-(t-Ts)/Td) * uu;
%         u(2) = (t-Ts)/Td * uu;
%     else
%         disp(["Error: LIP_func() t = " num2str(t)]);
%     end
%     
%     
%     dq(1,1) = q(2);
%     dq(2,1) = g/zc*q(1) + 1/m/zc*u(1) + 1/m/zc*u(2);
% end