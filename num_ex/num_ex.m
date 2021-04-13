clc; clear; close all;

%% Switched system eq.(44-45)
f0 = @(x) x - 1/(1+0.1*x^2);
f1 = @(x) x + 1/(1+0.1*x^2);

dt = 0.1;
t = 0:dt:100;
x = nan(length(t), 1);
y = nan(length(t), 1);

x(1) = 0;
y(1) = 0;

for i = 1:(length(t)-1)
    if y(i) == 0
        f = f0;
    else
        f = f1;
    end
    x(i+1) = f(x(i));
    y(i+1) = g(x(i), y(i));
    
end

figure(1);
plot(t, x, 'b', t, y, 'r');

Data.x = [x(1:end-1, :), y(1:end-1, :)];
Data.y = [x(2:end, :), y(2:end, :)];
Data.dt = 0.1;
save('data_ex.mat', 'Data');

%% Functions
% switched function
function lamb_kp1 = g(x_k, lamb_k)
    if x_k <= -1
        lamb_kp1 = 1;
    elseif x_k >= 1
        lamb_kp1 = 0;
    else
        lamb_kp1 = lamb_k;
    end
end

