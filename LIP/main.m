clc; 
clear;
close all;

%% Learn Koopman operator
if 0
tic;
load('data_mutistep.mat');
lamb = Data.x(:,5);

% Lift data to Koopman space
[psi_x_x_0, psi_x_y_0, psi_xu_x_0, psi_xu_y_0] = lift_data(Data.x(lamb==0, 1:4), Data.y(lamb==0, 1:4));
[psi_x_x_1, psi_x_y_1, psi_xu_x_1, psi_xu_y_1] = lift_data(Data.x(lamb==1, 1:4), Data.y(lamb==1, 1:4));

% Koopman operator (Eq. 17)
K0 = lsqminnorm([psi_x_x_0, psi_xu_x_0], psi_x_y_0);
K1 = lsqminnorm([psi_x_x_1, psi_xu_x_1], psi_x_y_1);

Dx = size(psi_x_x_0, 2);
Dxu = size(psi_xu_x_0, 2);

Kx0 = K0(1:Dx, :);
Kxu0 = K0(Dx+1:end, :);
Kx1 = K1(1:Dx, :);
Kxu1 = K1(Dx+1:end, :);

% Calculate A (Eq. 18)
% A0 = 1/Data.dt * logm(K0);
% A1 = 1/Data.dt * logm(K1);

toc;
save('Koopman.mat', 'Kx0', 'Kxu0', 'Kx1', 'Kxu1');

else
    load('data_mutistep.mat');
    load('Koopman.mat');
end

%% Estimate multistep
test = 1:size(Data.x_test,1);
x_test = Data.x_test;

x_est = estimate(Kx0, Kxu0, Kx1, Kxu1, x_test, Data.dt);

figure(1);
y_labels = {'x', 'Vx'};
% y_lims = [-0.4, 3; 0.75, 1; 1, 1.5; -1.5, 1.5];
for i = 1:2
    subplot(2,1,i); hold on;
    plot(test(x_test(test,5)==1)*Data.dt, x_test(x_test(test,5)==1,i), 'r.',...
        test(x_test(test,5)==0)*Data.dt, x_test(x_test(test,5)==0,i), 'b.');
    plot(test(x_est(test,5)==1)*Data.dt, x_est(x_est(test,5)==1,i), 'mx',...
        test(x_est(test,5)==0)*Data.dt, x_est(x_est(test,5)==0,i), 'cx', 'MarkerSize', 4);
    xlabel('t'); ylabel(y_labels{i});
%     ylim(y_lims(i,:));
end
legend('real single stance','real double stance','Koopman single stance','Koopman double stance');

figure(2);
xc_est = slit2continue(x_est);
xc_test = slit2continue(x_test);
error = xc_est(:, 1:2) - xc_test(:, 1:2);
rmse = rms(error, 1);
nrmse = rmse ./ (max(xc_est(:, 1:2), [], 1) - min(xc_est(:, 1:2), [], 1));
labels = categorical({'error x', 'error Vx'});
labels = reordercats(labels, {'error x', 'error Vx'});
bar(labels, nrmse);
