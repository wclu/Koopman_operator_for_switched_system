function out1 = func_M(t,in2)
%FUNC_M
%    OUT1 = FUNC_M(T,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    12-Jan-2021 11:10:07

Q1 = in2(1,:);
out1 = reshape([1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,6.0e1,0.0,0.0,0.0,0.0,Q1.^2.*6.0e1],[4,4]);
