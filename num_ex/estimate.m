% function x_est = estimate(K1, K0, g, x)
% function x_est = estimate(F1, F0, x)
function x_est = estimate(A1, A0, x)
% assume g is know now

x_est = x(1,:); % initial state
current_lamb = x(1,2);
current_state = x(1,1);
% f1 = @(t, x) f(F1, x);
% f0 = @(t, x) f(F0, x);
f1 = @(t, x) f(A1, x);
f0 = @(t, x) f(A0, x);
t_step = [0, 0.1]; % dt for the map
for i = 2:size(x, 1)
    if current_lamb ==1
        [~, psi] = ode45(f1, t_step, current_state);
    else
        [~, psi] = ode45(f0, t_step, current_state);
    end
    next_state = psi(end, 1); % Recover state from basis
    next_lamb = g(current_state, current_lamb);
    %next_lamb = x(i,2);
    x_est = cat(1, x_est, [next_state, next_lamb]);
    
    current_state = next_state;
    current_lamb = next_lamb;
end


end

% function dx = f(F, x)
%     % System representing Dubins' car.
%     [psi, ~, ~] = lift(x, x);
%     dx = F * psi;
% end

function dx = f(A, x)
    % System representing Dubins' car.
    dpsi = @(x) [0 1 2*x(1)]; % @(x)
    F = dpsi(x)' \ A';
    [psi, ~, ~] = lift(x, x);
    dx = F * psi;
end

%% switched function
function lamb_kp1 = g(x_k, lamb_k)
    if x_k >= 1
        lamb_kp1 = 1;
    elseif x_k <= 0
        lamb_kp1 = 0;
    else
        lamb_kp1 = lamb_k;
    end
end