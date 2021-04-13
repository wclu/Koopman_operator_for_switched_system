function [t, q, u, isdouble] = LIP_onestep(q0, p)
% p = {g, m, zc, q0, Ts, Td};
[g, m, zc, x0, Ts, Td] = deal(p{:});
kp = 300; kd = 100;
dt = 0.01;
[t_s, q_s] = ode45(@(t, q) LIP_single(t, q, p), 0:dt:(Ts-dt), q0);
[t_d, q_d] = ode45(@(t, q) LIP_double(t, q, p), Ts:dt:(Ts+Td-dt), q_s(end,:)');

t = [t_s; t_d(2:end,:)];
q = [q_s; q_d(2:end,:)];
isdouble = [false(length(t_s), 1); true(length(t_d)-1, 1)];

% Calculate input u
u = nan(length(t), 2);
for i = 1:length(t)
    traj = gen_traj(x0, Ts, Td, t(i));
    if ~isdouble(i,1)
        u(i,1) = [kp, kd] * (traj - q(i,:)');
        u(i,2) = 0;
    else
        uu = [kp, kd] * (traj - q(i,:)');
        u(i,1) = (1-(t(i)-Ts)/Td) * uu;
        u(i,2) = (t(i)-Ts)/Td * uu;
    end

end

end

%%
function dq = LIP_single(t, q, p)
    % q = [x, x_dot]'
    % u = [u1, u2]'
    % p = {g, m, zc, x0, Ts, Td};

    [g, m, zc, x0, Ts, Td] = deal(p{:});
    kp = 300; kd = 100;
    
    traj = gen_traj(x0, Ts, Td, t);
    u(1) = [kp, kd] * (traj - q);
    u(2) = 0;
    
    dq(1,1) = q(2);
    dq(2,1) = g/zc*q(1) + 1/m/zc*u(1) + 1/m/zc*u(2);
end

function dq = LIP_double(t, q, p)
    % q = [x, x_dot]'
    % u = [u1, u2]'
    % p = {g, m, zc, q0, Ts, Td};

    [g, m, zc, x0, Ts, Td] = deal(p{:});
    kp = 300; kd = 100;
    
    traj = gen_traj(x0, Ts, Td, t);
    uu = [kp, kd] * (traj - q);
    u(1) = (1-(t-Ts)/Td) * uu;
    u(2) = (t-Ts)/Td * uu;
    
    dq(1,1) = q(2);
    dq(2,1) = g/zc*q(1) + 1/m/zc*u(1) + 1/m/zc*u(2);
end