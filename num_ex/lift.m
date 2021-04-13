function [psi, psi_next, dpsi] = lift(x, y)
    % Function to lift a given data point x into Koopman space.
    % Basis for Koopman space shown below: (w=2)
    %
    % basis = [1 x1, x1^2]
    %
    % In addition, we need the derivative of the lifted data used in Eq 21
    % to calculate the nonlinear vector field. The Jacobian for the basis
    % is calculated analytically ans shown below:
    %
    % dpsi = [0 1 2*x(1)];

    % Set up basis coefficients for monomials
    [X1] = ndgrid(0:2);%%
    basis = [X1(:)];
    basis(sum(basis,2)>2,:) = [];
    %basis(sum(basis,2) == 0,:) = [];
    
    % Get psis to calculate Koopman Operator (Eq. 16)
    psi = x(1).^basis(:,1);
    psi_next = y(1).^basis(:,1);
    
    % Get analytical Jacobian for psi used in Eq. 21 
    dpsi = [0 1 2*x(1)]; % 2*x(1)
        
end