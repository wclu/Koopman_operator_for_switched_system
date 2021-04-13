clc; clear;

%% Parameters
g = 9.81;
l0 = 1;  
m = 60;  
k = 5600;
save('model/param.mat', 'm', 'g', 'l0', 'k');

%% Stance phase
syms t real
Q = sym('Q', [4, 1], 'real'); % Q=[q q_dot]'
q = Q(1:2); %q=[l phi]'
q_dot = Q(3:4); 

p = q(1) * [-sin(q(2)); cos(q(2))]; %p=[x y]'
J = jacobian(p, q);
p_dot = J * q_dot;

L = m * (p_dot' * p_dot) / 2 - (k * (q(1) - l0)^2 / 2 + m * g * p(2));

M = blkdiag(eye(2), jacobian(jacobian(L, q_dot), q_dot));
f = [q_dot; jacobian(L, q)' - jacobian(jacobian(L, q_dot), q) * q_dot];

matlabFunction(J, 'file', 'model/func_J', 'vars', {q});
matlabFunction(simplify(M), 'file', 'model/func_M', 'vars', {t, Q});
matlabFunction(f, 'file', 'model/func_f', 'vars', {t, Q});
matlabFunction(q(1) - l0, 1, 1, 'file', 'model/event_LO', ...
    'vars', {t, Q}, 'outputs', {'value', 'isterminal', 'direction'});

%%  Flight phase
syms beta real
P = sym('P', [4, 1], 'real');
p = P(1:2); p_dot = P(3:4);

h = [p_dot; 0; -g];
matlabFunction(h, 'file', 'model/func_h', 'vars', {t, P});
matlabFunction(p(2) - l0 * sin(beta), 1, -1, 'file', 'model/event_TD', ...
    'vars', {t, P, beta}, 'outputs', {'value', 'isterminal', 'direction'});

disp('Model done.');