function out1 = func_M(t,in2)
%FUNC_M
%    OUT1 = FUNC_M(T,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    16-Feb-2021 14:46:43

Q1 = in2(1,:);
out1 = reshape([1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,Q1.^2],[4,4]);
