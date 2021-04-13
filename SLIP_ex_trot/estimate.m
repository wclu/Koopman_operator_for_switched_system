function x_est = estimate(A1, A0, x, dt)
% assume g is know now

x_est = x(1,:); % initial state
current_lamb = x(1, 5);
current_state = x(1, 1:4);

% assume initial state is TD point
l0 = sqrt(x(1,1)^2 + x(1,2)^2);
phi0 = atan2(-x(1,1), x(1,2));

f1 = @(t, x) f(A1, x);
f0 = @(t, x) f(A0, x);
t_step = [0, dt]; % dt for the map
for i = 2:size(x, 1)
    if current_lamb ==1
        [~, xx] = ode45(f1, t_step, current_state);
    else
        [~, xx] = ode45(f0, t_step, current_state);
    end
    next_state = xx(end,:); % Recover state from basis
    [next_lamb, reset] = g(current_state, current_lamb, l0, phi0);
    %next_lamb = g(current_state, current_lamb, l0, phi0);
%     next_lamb = x(i,5);
    
    if reset
%         phi0 = acos(current_state(2)/l0);
%         current_state(1) = -sin(phi0);
        next_state(1) = -sin(phi0);
    end
    x_est = cat(1, x_est, [next_state, next_lamb]);
    
    current_state = next_state;
    current_lamb = next_lamb;
end


end

%% For test
% function x_est = estimate(A1, x, dt)
% % assume g is know now
% 
% x_est = x(1,1:4); % initial state
% % current_lamb = x(1, 5);
% current_state = x(1, 1:4);
% % f1 = @(t, x) f(F1, x);
% % f0 = @(t, x) f(F0, x);
% f1 = @(t, x) f(A1, x);
% % f0 = @(t, x) f(A0, x);
% t_step = [0, dt]; % dt for the map
% for i = 2:size(x, 1)
%     [~, xx] = ode45(f1, t_step, current_state);
%     next_state = xx(end,:); % Recover state from basis func_psi2x(psi);
%     %next_lamb = g(current_state, current_lamb);
%     x_est = cat(1, x_est, next_state);
%     
%     current_state = next_state;
% end
% 
% 
% end
%%
function dx = f(A, x)
    % System representing Dubins' car.
%     dpsi = @(x) [0 1 2*x(1) 0 x(2) 0 0 x(3) 0 0 0 x(4) 0 0 0;
%             0 0 0 1 x(1) 2*x(2) 0 0 x(3) 0 0 0 x(4) 0 0;
%             0 0 0 0 0 0 1 x(1) x(2) 2*x(3) 0 0 0 x(4) 0;
%             0 0 0 0 0 0 0 0 0 0 1 x(1) x(2) x(3) 2*x(4)];
%     F = dpsi(x)' \ A';
    F = func_dpsi(x)' \ A';
    [psi, ~] = lift(x, x);
    dx = F * psi;
end

% function y_kp1 = g(x_k, y_k, l0, phi0)
%     if (y_k == 1) && (sqrt(x_k(1)^2 + x_k(2)^2) > l0)
%         y_kp1 = 0; % LO
%     elseif (y_k == 0) && (x_k(2) < l0*cos(phi0))
%         y_kp1 = 1; % TD
%     else
%         y_kp1 = y_k;
%     end
% end

function [y_kp1, reset] = g(x_k, y_k, l0, phi0)
    if (y_k == 1) && (sqrt(x_k(1)^2 + x_k(2)^2) > l0)
        y_kp1 = 0; % LO
        reset = 0;
    elseif (y_k == 0) && (x_k(2) < l0*cos(phi0))
        y_kp1 = 1; % TD
        reset = 1;
    else
        y_kp1 = y_k;
        reset = 0;
    end
end