function dpsi_xu = func_dpsi_xu(in1,in2)
%FUNC_DPSI_XU
%    DPSI_XU = FUNC_DPSI_XU(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    22-Mar-2021 15:03:27

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1.^2;
t3 = x1.*x2.*2.0;
t4 = x2.^2;
t5 = u1.^2;
t6 = u1.*u2;
t7 = u2.^2;
t8 = u1.*x1.*2.0;
t9 = u1.*x2.*2.0;
t10 = u2.*x1;
t11 = u2.*x2;
t12 = x1.*x2;
t13 = u1.*x1;
t14 = u1.*x2;
t15 = u2.*x1.*2.0;
t16 = u2.*x2.*2.0;
t17 = u1.*u2.*2.0;
dpsi_xu = reshape([0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,x1.*2.0,0.0,0.0,0.0,t2.*3.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,x2,x1,0.0,0.0,t3,t2,0.0,0.0,0.0,x2.*2.0,0.0,0.0,t4,t3,0.0,0.0,0.0,t4.*3.0,0.0,0.0,0.0,0.0,1.0,0.0,u1,0.0,x1,0.0,t8,0.0,t2,0.0,0.0,u1,x2,0.0,t14,t13,t12,0.0,0.0,t9,t4,0.0,0.0,0.0,u1.*2.0,0.0,t5,0.0,t8,0.0,0.0,t5,t9,0.0,0.0,0.0,t5.*3.0,0.0,0.0,0.0,0.0,1.0,u2,0.0,0.0,x1,t15,0.0,0.0,t2,0.0,u2,0.0,x2,t11,t10,0.0,t12,0.0,t16,0.0,t4,0.0,0.0,u2,u1,t6,0.0,t10,t13,0.0,t6,t11,t14,0.0,0.0,t17,t5,0.0,0.0,0.0,u2.*2.0,t7,0.0,0.0,t15,0.0,t7,0.0,t16,0.0,0.0,t7,t17,0.0,0.0,0.0,t7.*3.0],[4,35]);
