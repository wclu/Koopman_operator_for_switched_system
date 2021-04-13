%% Fixed point test
clc; clear; close all;
addpath('model');

scale = pi / 180;
alpha_n = linspace(0, 30*scale, 121);
v = 1.25:0.25:2;
beta_lim = [55, 80];

%%
if 0
tic
for i = 1:length(v)
    beta_fxpnt = linspace(beta_lim(1, 1), beta_lim(1, 2), 121)*scale;%linspace(64, 71, 121)*scale;
    alpha_fxpnt = NaN(1,length(beta_fxpnt));
    flag_fxpnt = NaN(1,length(beta_fxpnt));

    for j = 1:length(beta_fxpnt)
        f = @(x) (sim_returnmap(sim_returnmap(x, v(i), beta_fxpnt(j)), v(i), beta_fxpnt(j)) - x);
        y = f(alpha_n); %alpha_np1-alpha_n;
        alpha_np1 = y+alpha_n;

        % returnmap ==============
        figure(i+1); hold on;
        plot(alpha_n/scale,alpha_np1/scale, 'marker','.','markersize',6);
        hold off;
        % ========================
        r = 1;
        for k = 1:1:length(y)-1
            if y(k)*y(k+1)<0 
                x1 = alpha_n(k);
                x2 = alpha_n(k+1);
                alpha_fxpnt(r,j) = fzero(f, [x1, x2]);
                if abs((alpha_np1(k)-alpha_np1(k+1))/(alpha_n(k)-alpha_n(k+1)))<1
                    flag_fxpnt(r,j) = 1;
                else
                    flag_fxpnt(r,j) = 0;
                end
                r = r+1;
            end
        end
        disp([i, j]);
    end
%     fxpnt = [beta_fxpnt; alpha_fxpnt]';
    save(['model\fxpnt_v', num2str(v(i) * 1000, '%d'), '_distribution.mat'], 'beta_fxpnt', 'alpha_fxpnt', 'flag_fxpnt')


    % returnmap ==============
    figure(i+1); hold on;
    plot([0 90],[0 90],'k');
    hold off; grid on;
    xlabel('\alpha_n (deg)');ylabel('\alpha_{n+1} (deg)');
    axis([0 90 0 90]);
    % ========================
end
toc
end

%%
figure(1)
hold on;
colormap jet;
map = colormap;
for i = 1:length(v)
    %===========
    load(['model\fxpnt_v', num2str(v(i) * 1000, '%d'), '_distribution.mat']);
    alpha_fxpnt = alpha_fxpnt(1,:);
    flag_fxpnt = flag_fxpnt(1,:);
%     alpha_fxpnt = reshape(alpha_fxpnt',1,[]);
%     flag_fxpnt = reshape(flag_fxpnt',1,[]);
%     beta_fxpnt = [beta_fxpnt, beta_fxpnt, beta_fxpnt, beta_fxpnt, beta_fxpnt];
%     RSLIP_beta_fxpnt = [RSLIP_beta_fxpnt, RSLIP_beta_fxpnt, RSLIP_beta_fxpnt, RSLIP_beta_fxpnt, RSLIP_beta_fxpnt];
    %===========
    n = int16(63/length(v)*i);
%     plot(beta_fxpnt/scale, alpha_fxpnt/scale, ...
%         'color',map(n,:), 'marker','o','markersize',3, 'MarkerFaceColor', map(n,:));
    plot(beta_fxpnt(flag_fxpnt==1)/scale, alpha_fxpnt(flag_fxpnt==1)/scale, ...
        'color',map(n,:), 'marker','o','markersize',3, 'MarkerFaceColor', map(n,:));
    plot(beta_fxpnt(flag_fxpnt==0)/scale, alpha_fxpnt(flag_fxpnt==0)/scale, ...
        'color',map(n,:), 'marker','o','markersize',3);
%     for j = 1:length(beta_fxpnt)   
%         if flag_fxpnt(j)==1
%             plot(beta_fxpnt(j)/scale, alpha_fxpnt(j)/scale, 'color',map(n,:), 'marker','o','markersize',3, 'MarkerFaceColor', map(n,:));
%         else
%             plot(beta_fxpnt(j)/scale, alpha_fxpnt(j)/scale, 'color',map(n,:), 'marker','o','markersize',3);
%         end
%         
%     end
end

xlabel('\beta_{SLIP} (deg)'); ylabel('\alpha (deg)');
legend(['v=' num2str(v(1))], ['v=' num2str(v(2))], ['v=' num2str(v(3))], ['v=' num2str(v(4))])
grid on;

%%
function alpha_np1 = sim_returnmap(alpha_n, v, beta)
    alpha_np1 = NaN(size(alpha_n));

    parfor k = 1:length(alpha_n)
        p_dot_TD = v * [cos(alpha_n(k)); -sin(alpha_n(k))];
        [~, ~, y, isstance] = SLIP_onestep(p_dot_TD, beta);
        alpha_np1(k) = atan2(-y(end, 4), y(end, 3));

        if alpha_np1(k) >= pi / 2 || all(isstance)
            alpha_np1(k) = NaN;
        end
    end
end