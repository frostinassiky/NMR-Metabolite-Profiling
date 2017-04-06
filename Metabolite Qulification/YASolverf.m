function [score,MINLP_output] = YASolverf(y,peak_part_th,peak_part_s,gauss_p1, pks_th,locs_th,locs_s,gauss_p2, pks_s, pstrength,peak_index)
meta_part = unique(peak_index(peak_part_th));
meta_order = zeros(max(meta_part),1);
for k = 1:length(meta_part)
    meta_order(meta_part(k)) = k;
end
% this should be var in future% Define variables
x = binvar(length(peak_part_th),length(peak_part_s));
% Define constraints 
[row,col,v] = find(pstrength==0);
for k = 1:length(v)
%     Constraints = [Constraints , x(row(k),col(k))==0];    
    x(row(k),col(k))=0;
end
Constraints = 0;
for k = 1:length(peak_part_th)
    Constraints = [Constraints , sum(x(k,:))<=1];    
end
for k = 1:length(peak_part_s)
    Constraints = [Constraints , sum(x(:,k))<=1];    
end
% and objective
for k = 1:length(peak_part_th)
for l = 1:length(peak_part_s)
	if pstrength(k,l) == 0
        continue;
    else
        pstrength(k,l) = gaussmf(locs_th(peak_part_th(k)),[gauss_p1, locs_s(peak_part_s(l))]) * ...
            gaussmf(pks_th(peak_part_th(k))*y(meta_order(peak_index(peak_part_th(k)))),[gauss_p2, pks_s(peak_part_s(l))]);
    end
end
end
Objective = -sum(x(:).*pstrength(:));
% Set some options for YALMIP and solver
ops = sdpsettings('verbose',0,'solver','bnb','showprogress',0);
% Solve the problem
sol = optimize(Constraints,Objective,ops);
if sol.problem == 0
    % Extract and display value
	[row,col] = find(value(x)==1);
    MINLP_output = struct('i',row','j',col','m',meta_part,'c',value(y));
else
    display('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
end
score = -value(Objective);
