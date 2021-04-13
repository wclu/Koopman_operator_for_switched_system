function [t, Q, P, isstance] = SLIP_onestep(p_dot_TD, beta)
%   Input:
%       p_dot_TD {2 x 1, real}:
%           Velocity at touchdown in Cartesian coordinates.
%       beta {1 x 1, real}:
%           Landing angle of the leg.
%
%   Output:
%       t {N x 1, real}:
%           Timestamps.
%       Q {N x 4, real}:
%           Q(k, :)' == [q(k); q_dot(k)]
%       P {N x 4, real}:
%           P(k, :)' == [p(k); p_dot(k)]
%       isstance {N x 1, logical}:
%           (isstance(k) == true) iff (in stance phase at t(k))

load('model/param.mat');
dt = 0.01;%eps^0.2;   % sampling period
reltol = eps^0.6;   abstol = eps^0.8;   % tolerances for ode45

% beta = pi/2 - phi0;

%% Stance phase
q0 = [l0; pi/2 - beta];
Q0 = [q0; func_J(q0) \ p_dot_TD];
opt_stance = odeset('event', @event_LO, 'mass', @func_M, ...
    'reltol', reltol, 'abstol', abstol);

[t_stance, Q_stance] = ode45(@func_f, 0:dt:1, Q0, opt_stance);
N_stance = length(t_stance);

% Map [q; q_dot] to [p; p_dot]
P_stance = NaN(N_stance, 4);
for i = 1:N_stance
    q = Q_stance(i, 1:2)';  q_dot = Q_stance(i, 3:4)';
    P_stance(i, :) = [q(1) * [-sin(q(2)); cos(q(2))]; func_J(q) * q_dot]';
end

%% Flight phase
delta_p2 = l0 * sin(beta) - P_stance(end, 2);
p2_dot_LO = P_stance(end, 4);
% Estimate flight duration
T = (p2_dot_LO + sqrt(p2_dot_LO^2 - 2 * g * delta_p2)) / g;

if isreal(T) && (T >= dt)
    t_LO = t_stance(end); y0 = P_stance(end, :)';
    opt_flight = odeset('event', @(t, y) event_TD(t, y, beta), ...
        'reltol', reltol, 'abstol', abstol);

    [t_flight, y_flight] = ode45(@func_h, t_LO + (0:dt:1), y0, opt_flight);
    
    % if eventTD not work ===
     if t_flight(end) >= t_LO+1-dt
        for j = length(t_flight):-1:1
            if (l0*sin(beta)-y_flight(j,2))<0
                t_flight = t_flight(1:j);
                y_flight = y_flight(1:j,:);
                break;
            end
        end
     end
    %========================
    N_flight = length(t_flight);
    phi0 = pi/2 - beta;
    Q_fight = [l0*ones(N_flight - 1, 1), phi0*ones(N_flight - 1, 1),...
        zeros(N_flight - 1, 2)];

    t = [t_stance; t_flight(2:end)];
    Q = [Q_stance; Q_fight];%NaN(N_flight - 1, 4)
    P = [P_stance; y_flight(2:end, :)];  
    isstance = [true(N_stance, 1); false(N_flight - 1, 1)];

else
    warning('!! Failed to liftoff !!')
    t = t_stance;
    Q = Q_stance;
    P = P_stance;
    isstance = true(N_stance, 1);
end

end