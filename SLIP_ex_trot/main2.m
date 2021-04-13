clc; 
clear;
close all;

%% Learn Koopman operator
load('data_mutistep.mat');
% lamb = Data.x(:,5);

% Lift data to Koopman space
[psi_x, psi_y] = lift_data(Data.x(:, 1:4), Data.y(:, 1:4));

% Koopman operator (Eq. 17)
K = lsqminnorm(psi_x, psi_y);

% Calculate A (Eq. 18)
A = 1/Data.dt * logm(K);

%% Estimate multistep
test = 1:size(Data.x_test1,1);
x_test = Data.x_test1;
% test = 1:size(Data.x_test2,1);
% x_test = Data.x_test2;
x_est = estimate2(A, x_test, Data.dt);

figure(1);
y_labels = {'x', 'y', 'Vx', 'Vy'};
% y_lims = [-0.4, 3; 0.75, 1; 1, 1.5; -1.5, 1.5];
% y_lims = [-0.4, 3; 0.75, 1; 1, 1.5; -1.5, 1.5];
for i = 1:4
    subplot(2,2,i); hold on;
    plot(test(x_test(test,5)==1)*Data.dt, x_test(x_test(test,5)==1,i), 'r.',...
        test(x_test(test,5)==0)*Data.dt, x_test(x_test(test,5)==0,i), 'b.');
    plot(test*Data.dt, x_est(:,i), 'mx', 'MarkerSize', 4);
%     plot(test(x_est(test,5)==1)*Data.dt, x_est(x_est(test,5)==1,i), 'mx',...
%         test(x_est(test,5)==0)*Data.dt, x_est(x_est(test,5)==0,i), 'cx', 'MarkerSize', 4);
    xlabel('t'); ylabel(y_labels{i});
%     ylim(y_lims(i,:));
end
legend('real stance','real flight','Koopman');
%%
figure(2);
xc_est = slit2continue(x_est);
xc_test = slit2continue(x_test);
error = xc_est(:, 1:4) - xc_test(:, 1:4);
rmse = rms(error, 1);
nrmse = rmse ./ (max(xc_est(:, 1:4), [], 1) - min(xc_est(:, 1:4), [], 1));
labels = categorical({'error x', 'error y', 'error Vx', 'error Vy'});
labels = reordercats(labels, {'error x', 'error y', 'error Vx', 'error Vy'});
bar(labels, nrmse);