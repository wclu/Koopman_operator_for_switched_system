function q = gen_traj(q0, Ts, Td, t)
% Input:
%   Ts: duration of single stance
%   Td: duration of double stance
%   q0: initial condition [x0, x_dot0]
%   t: time (scalar)
% Return:
%   q: state [x, x_dot] at time t

    % clc; clear; close all;
    % Ts = 1;
    % Td = 0.5;
    % q0 = [-0.5, 2]';
    
    t = mod(t, Ts+Td);

    if t>=0 && t<Ts
        % single stance
        A = [0, 0, 0, 1;
             Ts^3, Ts^2, Ts, 1;
             0, 0, 1, 0;
             3*Ts^2, 2*Ts, 1, 0];
        b = [q0(1); -q0(1); q0(2); q0(2)];
        a = A\b;
        q = [polyval(a, t);
             polyval(polyder(a), t)];

    elseif t>=Ts && t<=(Ts+Td)
        % double stance
        q = [polyval([q0(2), -q0(1)], t-Ts);
        	 ones(size(t)) * q0(2)];
    else
        disp(["Error: gen_traj() t = " num2str(t)]);
    end

end
%% Plot and check
% t = 0:0.01:Ts;
% x = polyval(a, t);
% dx = polyval(polyder(a), t);
% ddx = polyval(polyder(polyder(a)), t);
% figure(1);
% subplot(3,1,1), plot(t, x, 'b'); hold on;
% subplot(3,1,2), plot(t, dx, 'b'); hold on;
% subplot(3,1,3), plot(t, ddx, 'b'); hold on;
% 
% xf = x(end);
% dxf = dx(end);
% 
% t = Ts:0.01:(Ts+Td);
% x = polyval([dxf, xf], t-Ts);
% dx = ones(size(t)) * dxf;
% ddx = zeros(size(t));
% figure(1);
% subplot(3,1,1), plot(t, x, 'r');
% subplot(3,1,2), plot(t, dx, 'r');
% subplot(3,1,3), plot(t, ddx, 'r');