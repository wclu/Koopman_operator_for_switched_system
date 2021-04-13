clc; 
clear;
close all;

%%
% load('data_ex.mat');
load('data_ex2.mat');
lamb = Data.x(:,2);


% Lift data to Koopman space
[psi_x_1, psi_y_1, dpsi_1] = lift_data(Data.x(lamb==1, 1), Data.y(lamb==1, 1));
[psi_x_0, psi_y_0, dpsi_0] = lift_data(Data.x(lamb==0, 1), Data.y(lamb==0, 1));

% Koopman operator (Eq. 17)
% K1 = pinv(psi_x_1) * psi_y_1;
% K0 = pinv(psi_x_0) * psi_y_0;
K1 = lsqminnorm(psi_x_1, psi_y_1);
K0 = lsqminnorm(psi_x_0, psi_y_0);

% Calculate A (Eq. 18)
A1 = 1/Data.dt * logm(K1);
A0 = 1/Data.dt * logm(K0);


% Estimate
test = 1:500;
x_test = Data.x(test, :);
x_est = estimate(A1, A0, x_test);

figure(1); hold on;
plot(test, x_test(:,1), 'b', test, x_test(:,2), 'k');
plot(test, x_est(:,1), 'r--', test, x_est(:,2), 'm--');%
hold off;