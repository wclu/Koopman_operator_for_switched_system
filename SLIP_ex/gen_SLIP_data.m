clc; clear; close all;
addpath('model');

%% Setting
alpha_list = (10:5:20)/180*pi;
v_list = 1.5:0.25:2;

beta_lim = [63, 70; 59, 68; 56, 67;
    64, 71; 61, 69; 58, 68;
    66, 72; 63, 71;60, 69]/180*pi;
beta_list = nan(9,5);
for i = 1:size(beta_list, 1)
    beta_list(i,:) = linspace(beta_lim(i, 1), beta_lim(i, 2), 5);
end

% Data.t = [];
% Data.Q = [];
% Data.P = [];
% Data.isstance = [];
Data.x = []; % input x_i
Data.y = []; % output x_i+1
Data.x_test = [];
Data.dt = 0.01;

%% Generate training data
for i = 1:length(alpha_list)
    for j = 1:length(v_list)
        for k = 1:size(beta_list, 2)
            alpha = alpha_list(i);
            v = v_list(j);
            beta = beta_list((i-1)*length(v_list)+j, k);
            p_dot_TD = v*[cos(alpha); -sin(alpha)];
            
            [t, Q, P, isstance] = SLIP_onestep(p_dot_TD, beta);
            %Data.t = [Data.t; t];
            %Data.Q = [Data.Q; Q];
            %Data.P = [Data.P; P];
            %Data.isstance = [Data.isstance; isstance];
            %P(:,1) = P(:,1) - P(1,1);
            Data.x = [Data.x; [P(1:end-1, :), isstance(1:end-1, 1)]];
            Data.y = [Data.y; [P(2:end, :), isstance(2:end, 1)]];
            
        end
    end
end

% save('data.mat', 'Data');

%% test
alpha = 0.2117; %0.2060, 0.2117
p_dot_TD = 1.5*[cos(alpha); -sin(alpha)];
beta = 1.2392; %1.2382, 1.2392

% x0 = 0;
for i = 1:3
    [t, Q, P, isstance] = SLIP_onestep(p_dot_TD, beta);

    figure(1);
    plot(P(isstance,1), P(isstance,2), 'r', ...
        P(~isstance,1), P(~isstance,2), 'b');
    hold on;
    
    p_dot_TD = P(end, 3:4)';
%     P = P + [x0, 0, 0, 0];
%     x0 = P(end, 1);
    Data.x_test = [Data.x_test; [P(1:end-1, :), isstance(1:end-1, 1)]];
end
Data.x_test = [Data.x_test; [P(end, :), isstance(end, 1)]];

save('data.mat', 'Data');

%% Fixed point test
% % result: (v, alpha, beta) = (1.5, 0.2060, 1.2382), (1.5, 0.2117, 1.2392)
% scale = pi / 180;
% alpha_n = linspace(0, 20*scale, 121);
% v = 1.5;
% beta = linspace(64, 71, 121)*scale;
% alpha_fxpnt = NaN(5,length(beta));
% flag_fxpnt = NaN(5,length(beta));
% 
% for j = 1:length(beta)
%     f = @(x) (sim_returnmap(sim_returnmap(x, v, beta(j)), v, beta(j)) - x);
%     y = f(alpha_n); %alpha_np1-alpha_n;
%     alpha_np1 = y+alpha_n;
% 
%     % returnmap ==============
%     figure(1); hold on;
%     plot(alpha_n/scale,alpha_np1/scale, 'marker','.','markersize',6);
%     hold off;
%     % ========================
%     r = 1;
%     for k = 1:1:length(y)-1
%         if y(k)*y(k+1)<0 
%             x1 = alpha_n(k);
%             x2 = alpha_n(k+1);
%             alpha_fxpnt(r,j) = fzero(f, [x1, x2]);
%             if abs((alpha_np1(k)-alpha_np1(k+1))/(alpha_n(k)-alpha_n(k+1)))<1
%                 flag_fxpnt(r,j) = 1;
%             else
%                 flag_fxpnt(r,j) = 0;
%             end
%             r = r+1;
%         end
%     end
% end
% 
% % returnmap ==============
% figure(1); hold on;
% plot([0 90],[0 90],'k');
% hold off; grid on;
% xlabel('\alpha_n (deg)');ylabel('\alpha_{n+1} (deg)');
% axis([0 90 0 90]);
% % ========================
% 
% function alpha_np1 = sim_returnmap(alpha_n, v, beta)
%     alpha_np1 = NaN(size(alpha_n));
% 
%     parfor k = 1:length(alpha_n)
%         p_dot_TD = v * [cos(alpha_n(k)); -sin(alpha_n(k))];
%         [~, ~, y, isstance] = SLIP_onestep(p_dot_TD, beta);
%         alpha_np1(k) = atan2(-y(end, 4), y(end, 3));
% 
%         if alpha_np1(k) >= pi / 2 || all(isstance)
%             alpha_np1(k) = NaN;
%         end
%     end
% end