function [psi, psi_next] = lift(x, y)
% function [psi, psi_next, dpsi] = lift(x, y)
    % Function to lift a given data point x into Koopman space.
    % Basis for Koopman space shown below: (w=2)
    %
    % basis = [1 x1, x1^2, x2, x1x2, x2^2, x3, x1x3, x2x3, x3^2, x4, x1x4, x2x4, x3x4, x4^2]
    %
    % In addition, we need the derivative of the lifted data used in Eq 21
    % to calculate the nonlinear vector field. The Jacobian for the basis
    % is calculated analytically ans shown below:
    %
    % dpsi = [0 1 2*x(1) 0 x(2) 0 0 x(3) 0 0 0 x(4) 0 0 0;
    %         0 0 0 1 x(1) 2*x(2) 0 0 x(3) 0 0 0 x(4) 0 0;
    %         0 0 0 0 0 0 1 x(1) x(2) 2*x(3) 0 0 0 x(4) 0;
    %         0 0 0 0 0 0 0 0 0 0 1 x(1) x(2) x(3) 2*x(4)];
%=====
%     % Set up basis coefficients for monomials
%     [X1,X2,X3,X4] = ndgrid(0:2);
%     basis = [X1(:),X2(:),X3(:),X4(:)];
%     basis(sum(basis,2)>2,:) = [];
%     %basis(sum(basis,2) == 0,:) = [];
%     
%     % Get psis to calculate Koopman Operator (Eq. 16)
%     psi = x(1).^basis(:,1) .* x(2).^basis(:,2) .* x(3).^basis(:,3) .* x(4).^basis(:,4);
%     psi_next = y(1).^basis(:,1) .* y(2).^basis(:,2) .* y(3).^basis(:,3) .* y(4).^basis(:,4);
%     
%     % Get analytical Jacobian for psi used in Eq. 21 
% %     dpsi = [0 1 2*x(1) 0 x(2) 0 0 x(3) 0 0 0 x(4) 0 0 0;
% %             0 0 0 1 x(1) 2*x(2) 0 0 x(3) 0 0 0 x(4) 0 0;
% %             0 0 0 0 0 0 1 x(1) x(2) 2*x(3) 0 0 0 x(4) 0;
% %             0 0 0 0 0 0 0 0 0 0 1 x(1) x(2) x(3) 2*x(4)];
%=====
    % x, y are Nx1 or 1xN 
    % psi is Mx1
    if (size(x,1)==1) && (size(x,2)~=1)
        x = x';
        y = y';
    end
    
    psi = func_psi(x);
    psi_next = func_psi(y);
    
end