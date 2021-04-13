clc; clear;

%% Setting
x = sym('x', [2, 1], 'real');   % states
u = sym('u', [2, 1], 'real');   % system input
w = 3;                          % max nominal power

%% Generate nominal basis
% For psi_x(x_k)
n = size(x,1);    
basis = gen_basis(n,w);

psi_x = 1;
for i = 1:n
    psi_x = psi_x .* x(i) .^ basis(:,i);
end
dpsi_x = jacobian(psi_x, x)';

disp('psi_x.T = '); disp(psi_x');
disp('dpsi_x = '); disp(dpsi_x);

matlabFunction(psi_x, 'file', 'func_psi_x', 'vars', {x});
matlabFunction(dpsi_x, 'file', 'func_dpsi_x', 'vars', {x});

% For psi_xu(x_k, u_k)
if ~isnan(u)
    n = size(x,1) + size(u,1);
    xx = cat(1, x, u);
    
    basis = gen_basis(n,w);

    psi_xu = 1;
    for i = 1:n
        psi_xu = psi_xu .* xx(i) .^ basis(:,i);
    end
    dpsi_xu = jacobian(psi_xu, xx)';

    disp('psi_xu.T = '); disp(psi_xu');
    disp('dpsi_xu = '); disp(dpsi_xu);

    matlabFunction(psi_xu, 'file', 'func_psi_xu', 'vars', {x,u});
    matlabFunction(dpsi_xu, 'file', 'func_dpsi_xu', 'vars', {x,u});
end

%%
function basis = gen_basis(n,w)
    X = cell(1, n);
    [X{:}] = ndgrid(0:w);
    basis = [];
    for i = 1:n
        basis = cat(2, basis, reshape(X{i}, [], 1));
    end
    basis(sum(basis,2)>w,:) = [];
end