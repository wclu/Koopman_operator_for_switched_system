function dpsi_x = func_dpsi_x(in1)
%FUNC_DPSI_X
%    DPSI_X = FUNC_DPSI_X(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    22-Mar-2021 15:03:26

x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1.^2;
t3 = x1.*x2.*2.0;
t4 = x2.^2;
dpsi_x = reshape([0.0,0.0,1.0,0.0,x1.*2.0,0.0,t2.*3.0,0.0,0.0,1.0,x2,x1,t3,t2,0.0,x2.*2.0,t4,t3,0.0,t4.*3.0],[2,10]);