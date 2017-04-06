function [x,fval,exitflag,output,population,score] = optimal_ga(nvars,lb,ub,fun)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'StallGenLimit', 50000);
options = gaoptimset(options,'TolFun', 1e-12);
options = gaoptimset(options,'TolCon', 1e-6);
options = gaoptimset(options,'Display', 'iterative');
options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv });
[x,fval,exitflag,output,population,score] = ...
ga(@(x)norm(fun(x),2),nvars,[],[],[],[],lb,ub,[],[],options);
