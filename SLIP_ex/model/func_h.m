function h = func_h(t,in2)
%FUNC_H
%    H = FUNC_H(T,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    12-Jan-2021 11:10:07

P3 = in2(3,:);
P4 = in2(4,:);
h = [P3;P4;0.0;-9.81e2./1.0e2];