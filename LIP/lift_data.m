function [psi_x_x, psi_x_y, psi_xu_x, psi_xu_y] = lift_data(x, y)
    % Purpose of this function is to lift the entire sequence of data from
    % the simulation into Koopman space and return the psi vectors (Eq. 16)
    
    psi_x_x = [];                     % Input psi vector (Eq.16)
    psi_x_y = [];                     % Output psi vector (Eq.16)
    psi_xu_x = [];                    % Input psi vector (Eq.16)
    psi_xu_y = [];                    % Output psi vector (Eq.16)
    
    for i = 1:size(x,1)
        % Lift to Koopman space
        [psi_x, psi_x_next, psi_xu, psi_xu_next] = lift(x(i,:), y(i,:));
             
        psi_x_x = cat(1, psi_x_x, psi_x');            % Concatenate input psi
        psi_x_y = cat(1, psi_x_y, psi_x_next');       % Concatenate output psi
        psi_xu_x = cat(1, psi_xu_x, psi_xu');         % Concatenate input psi
        psi_xu_y = cat(1, psi_xu_y, psi_xu_next');    % Concatenate output psi
  
    end
    
    % For Koopman Operator calculation, we need the dimensions of the psi
    % matrix to agree with the other matrices used in Eq. 21, so we
    % truncate.
    samples1 = size(psi_x_x,1);                % Number of rows in psi
    dim1 = size(psi_x_x, 2);                   % Number of sim samples
    k1 = dim1*floor(samples1/dim1);            % Rescaling factor
    samples2 = size(psi_x_x,1);                % Number of rows in psi
    dim2 = size(psi_x_x, 2);                   % Number of sim samples
    k2 = dim2*floor(samples2/dim2);            % Rescaling factor
    
    psi_x_x = psi_x_x(1:k1,:);                 % Rescale input psi
    psi_x_y = psi_x_y(1:k1,:);                 % Rescale output psi
    psi_xu_x = psi_xu_x(1:k2,:);               % Rescale input psi
    psi_xu_y = psi_xu_y(1:k2,:);               % Rescale output psi
    
end