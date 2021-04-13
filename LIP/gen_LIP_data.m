clc; clear; close all;

%% Setting
% model parameter
g = 9.81;
m = 1;
zc = 1;

% gait parameter
q0_list = [-0.5, -0.4, -0.3];
dq0_list = [2.3, 2.5, 2.7];
Ts_list = [0.7, 0.8, 0.9];
Td_list = [0.3, 0.4, 0.5];

% test case
test_idx = [2, 2, 2, 2];
% p_test = {g, m, zc, [q0_list(2); dq0_list(2)], Ts_list(2), Td_list(2)};

% save data
Data.x = [];
Data.y = [];
Data.x_test = [];
Data.dt = 0.01;

%% Generate data
if 1
tic;
for i = 1:length(q0_list)
    for j = 1:length(dq0_list)
        for k = 1:length(Ts_list)
            for l = 1:length(Td_list)
                q0 = [q0_list(i); dq0_list(j)];
                x0 = 0; x = [];
                t0 = 0;
                for n = 1:2
                    
                    p = {g, m, zc, q0, Ts_list(k), Td_list(l)};
                    [t, q, u, isdouble] = LIP_onestep(q0, p);
                    
                    % reset q0, Note: q0(1,1) = q0(1,1)
                    q0(2,1) = q(end, 2);
                    
                    % continuous trajectory
                    if sum([i, j, k, l]==test_idx)~=length(test_idx)
                        if n==1
                            x0 = q(1,1);
                        end
                        q = q + [x0 - q(1,1), 0];
                        t = t + t0;
                        x0 = q(end, 1);
                        t0 = t(end);
                    end
                    
                    figure(1);
                    subplot(4,1,1); hold on;
                    plot(t(~isdouble), q(~isdouble,1), 'b', t(isdouble), q(isdouble,1), 'r');
                    xlabel('t'); ylabel('x');
                    subplot(4,1,2); hold on;
                    plot(t(~isdouble), q(~isdouble,2), 'b', t(isdouble), q(isdouble,2), 'r');
                    xlabel('t'); ylabel('Vx');
                    subplot(4,1,3); hold on;
                    plot(t(~isdouble), u(~isdouble,1), 'b', t(isdouble), u(isdouble,1), 'r');
                    xlabel('t'); ylabel('u1');
                    subplot(4,1,4); hold on;
                    plot(t(~isdouble), u(~isdouble,2), 'b', t(isdouble), u(isdouble,2), 'r');
                    xlabel('t'); ylabel('u2');
                    
                    x = [x; [q(1:end-1, :), u(1:end-1, :), isdouble(1:end-1, 1)]];
                end
                
                if sum([i, j, k, l]==test_idx)~=length(test_idx)
                    Data.x = [Data.x; x(1:end-1, :)];
                    Data.y = [Data.y; x(2:end, :)];
                else
                    Data.x_test = x(1:end-1, :);
                end
            end
        end
    end
end
save('data_mutistep.mat', 'Data');
toc;
end


