function x_est = estimate(Kx0, Kxu0, Kx1, Kxu1, X, dt)
% assume g is know now

x_est = X(1,:); % initial state
current_lamb = X(1, 5);
current_state = X(1, 1:2);
u = X(:, 3:4);

% Reset parameter
x0 = X(1,1);
% reset = 0;

% t_step = [0, dt]; % dt for the map
for i = 2:size(X, 1)
    current_lamb = X(i-1, 5);
    next_lamb = X(i, 5);
    if current_lamb ==0
        next_state = f(Kx0, Kxu0, current_state, u(i-1,:));
%         f0 = @(t, x) f(Kx0, Kxu0, x, u(i-1,:));
%         [~, xx] = ode45(f0, t_step, current_state);
    else
        next_state = f(Kx1, Kxu1, current_state, u(i-1,:));
%         f1 = @(t, x) f(Kx1, Kxu1, x, u(i-1,:));
%         [~, xx] = ode45(f1, t_step, current_state);
    end
%     next_state = xx(end,:);
%     [next_lamb, reset] = g(current_state, current_lamb, l0, phi0);
    %next_lamb = g(current_state, current_lamb, l0, phi0);
%     next_lamb = x(i,5);
    
    reset = (current_lamb==1) && (next_lamb==0);
    if  reset % reset
        next_state(1) = x0;
    end
    
    x_est = cat(1, x_est, [next_state, u(i,:), next_lamb]);
    current_state = next_state;
end


end

%%
function x_kp1 = f(Kx, Kxu, x_k, u)
    xx = [x_k, u];
    [psi_x, ~, psi_xu, ~] = lift(xx, xx);
    psi_y = (psi_x'*Kx + psi_xu'*Kxu)';
    x_kp1 = psi_y([2,5])'; % states from nominal
end
% function dx = f(A, x)
%     F = func_dpsi(x)' \ A';
%     [psi, ~] = lift(x, x);
%     dx = F * psi;
% end

% function y_kp1 = g(x_k, y_k, l0, phi0)
%     if (y_k == 1) && (sqrt(x_k(1)^2 + x_k(2)^2) > l0)
%         y_kp1 = 0; % LO
%     elseif (y_k == 0) && (x_k(2) < l0*cos(phi0))
%         y_kp1 = 1; % TD
%     else
%         y_kp1 = y_k;
%     end
% end

% function [y_kp1, reset] = g(x_k, y_k, l0, phi0)
%     if (y_k == 1) && (sqrt(x_k(1)^2 + x_k(2)^2) > l0)
%         y_kp1 = 0; % LO
%         reset = 0;
%     elseif (y_k == 0) && (x_k(2) < l0*cos(phi0))
%         y_kp1 = 1; % TD
%         reset = 1;
%     else
%         y_kp1 = y_k;
%         reset = 0;
%     end
% end