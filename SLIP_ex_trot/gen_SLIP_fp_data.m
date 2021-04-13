clc; clear; close all;
addpath('model');

%% Setting
alpha_list = (12:4:24)/180*pi;
v_list = 0.95:0.025:1.05;

test = [v_list(2), alpha_list(2); v_list(4), alpha_list(3)];

Data.x = []; % input x_i
Data.y = []; % output x_i+1
Data.x_test1 = [];
Data.x_test2 = [];
Data.dt = 0.01;

fp = [];
% x_dot, y_dot, alphaTD

%% Generate training & testing data
for i = 1:length(v_list)
    fxpnt = load_fxpnt(v_list(i));
    beta_fxpnt = fxpnt.beta_fxpnt;
    alpha_fxpnt = fxpnt.alpha_fxpnt(1,:);
    
    beta_list = spline(alpha_fxpnt, beta_fxpnt, alpha_list);
    
    for j = 1:length(alpha_list)
        p_dot_TD = v_list(i)*[cos(alpha_list(j)); -sin(alpha_list(j))];
        x0 = 0;
        x = [];
        
        % ===
        fp = [fp; [p_dot_TD', beta_list(j)]];
        % ===
        
        if v_list(i)==test(1,1) && alpha_list(j)==test(1,2)
            s = 3;
        elseif v_list(i)==test(2,1) && alpha_list(j)==test(2,2)
            s = 3;
        else
            s = 3;
        end
        
        for k = 1:s
            [t, Q, P, isstance] = SLIP_onestep(p_dot_TD, beta_list(j));
            
            p_dot_TD = P(end, 3:4)';
            
            if (v_list(i)~=test(1,1) || alpha_list(j)~=test(1,2)) &&...
               (v_list(i)~=test(2,1) || alpha_list(j)~=test(2,2))
                if k==1
                    x0 = P(1,1);
                end
                P = P + [x0 - P(1,1), 0, 0, 0];%
                x0 = P(end, 1);
%                 figure(1);
%                 plot(P(isstance,1), P(isstance,2), 'r', ...
%                     P(~isstance,1), P(~isstance,2), 'b');
%                 hold on;
%                 disp([v_list(i), alpha_list(j)*180/pi]);
            end

            figure(1);
            plot(P(isstance,1), P(isstance,2), 'r', ...
                P(~isstance,1), P(~isstance,2), 'b');
            hold on;

            x = [x; [P(1:end-1, :), isstance(1:end-1, 1)]];
        end
        
        % Save data points
        if v_list(i)==test(1,1) && alpha_list(j)==test(1,2)
            Data.x_test1 = x(1:end-1, :);
            %disp(['(v, alpha) = (', num2str(v_list(i)), ', ', num2str(alpha_list(j)), ')']);
        elseif v_list(i)==test(2,1) && alpha_list(j)==test(2,2)
            Data.x_test2 = x(1:end-1, :);
            %disp(['(v, alpha) = (', num2str(v_list(i)), ', ', num2str(alpha_list(j)), ')']);
        else
            Data.x = [Data.x; x(1:end-1, :)];
            Data.y = [Data.y; x(2:end, :)];
        end
        

    end
end

save('data_mutistep.mat', 'Data');
save('fp_trot.mat', 'fp');

%% 
function fxpnt = load_fxpnt(v)
    if any(1.25:0.25:2)
        fxpnt = load(['fxpnt_v', num2str(1000 * v), '_distribution.mat']);
    else
        error('Velocity can only be 1.25:0.25:2 [m/s]')
    end
end