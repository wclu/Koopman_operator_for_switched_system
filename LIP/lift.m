function [psi_x, psi_x_next, psi_xu, psi_xu_next] = lift(x, y)
    % x, y are Nx1 or 1xN 
    % psi is Mx1
    if (size(x,1)==1) && (size(x,2)~=1)
        x = x';
        y = y';
    end
    
    psi_x = func_psi_x(x(1:2));
    psi_x_next = func_psi_x(y(1:2));
    psi_xu = func_psi_xu(x(1:2), x(3:4));
    psi_xu_next = func_psi_xu(y(1:2), y(3:4));
    
end