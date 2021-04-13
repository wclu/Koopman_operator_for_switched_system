function psi_xu = func_psi_xu(in1,in2)
%FUNC_PSI_XU
%    PSI_XU = FUNC_PSI_XU(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    22-Mar-2021 15:03:26

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1.^2;
t3 = x2.^2;
t4 = u1.^2;
t5 = u2.^2;
psi_xu = [1.0;x1;t2;t2.*x1;x2;x1.*x2;t2.*x2;t3;t3.*x1;t3.*x2;u1;u1.*x1;t2.*u1;u1.*x2;u1.*x1.*x2;t3.*u1;t4;t4.*x1;t4.*x2;t4.*u1;u2;u2.*x1;t2.*u2;u2.*x2;u2.*x1.*x2;t3.*u2;u1.*u2;u1.*u2.*x1;u1.*u2.*x2;t4.*u2;t5;t5.*x1;t5.*x2;t5.*u1;t5.*u2];