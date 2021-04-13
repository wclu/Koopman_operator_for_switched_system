clc; clear;

%% Setting
x = sym('x', [4, 1], 'real');   % states
u = nan;                        % system input
w = 3;                          % max nominal power

%% Generate nominal basis
if isnan(u)
    n = size(x,1);
    xx = x;
else
    n = size(x,1) + size(u,1);
    xx = cat(1, x, u);
end

X = cell(1, n);
[X{:}] = ndgrid(0:w);
basis = [];
for i = 1:n
    basis = cat(2, basis, reshape(X{i}, [], 1));
end
basis(sum(basis,2)>w,:) = [];

psi = 1;
psi_next = 1;
for i = 1:n
    psi = psi .* xx(i) .^ basis(:,i);
end
% psi = psi';
% psi_next = psi_next';
dpsi = jacobian(psi, xx)';

disp('psi.T = '); disp(psi');
disp('psi_next.T = '); disp(psi_next');
disp('dpsi = '); disp(dpsi);

%% Save functions
if isnan(u)
    matlabFunction(psi, 'file', 'func_psi', 'vars', {x});
    matlabFunction(dpsi, 'file', 'func_dpsi', 'vars', {x});
else
    matlabFunction(psi, 'file', 'func_psi', 'vars', {x, u});
    matlabFunction(dpsi, 'file', 'func_dpsi', 'vars', {x, u});
end

% idx = find(sum(basis,2)==1);
% z = sym('z', size(psi), 'real');
% states = z(idx,1);
% matlabFunction(states, 'file', 'func_psi2x', 'vars', {z});