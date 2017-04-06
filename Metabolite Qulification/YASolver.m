meta_part = unique(peak_index(peak_part_th));
% this should be var in future% Define variables
x = binvar(length(peak_part_th),length(peak_part_s));
y = sdpvar(settings.num,1);
% Define constraints 
Constraints = y>=0;
[row,col,v] = find(pstrength==0);
for k = 1:settings.num
%     Constraints = [Constraints , x(row(k),col(k))==0];
    if k == meta_part
        continue
    end
    y(k) = 0;
end
for k = 1:length(v)
%     Constraints = [Constraints , x(row(k),col(k))==0];    
    x(row(k),col(k))=0;
end
for k = 1:length(peak_part_th)
    Constraints = [Constraints , sum(x(k,:))<=1];    
end
for k = 1:length(peak_part_s)
    Constraints = [Constraints , sum(x(:,k))<=1];    
end
% and objective
pstrength_var = sdpvar(length(peak_part_th),length(peak_part_s),'full');
for k = 1:length(peak_part_th)
for l = 1:length(peak_part_s)
	if pstrength(k,l) == 0
        pstrength_var(k,l) = 0;
    else
        pstrength_var(k,l) = gaussmf(locs_th(peak_part_th(k)),[gauss_p1, locs_s(peak_part_s(l))]) * ...
            gaussmf(pks_th(peak_part_th(k))*y(peak_index(peak_part_th(k))),[gauss_p2, pks_s(peak_part_s(l))]);
    end
end
end
Objective = -sum(x(:).*pstrength_var(:));
% Set some options for YALMIP and solver
ops = sdpsettings('verbose',1,'solver','bnb');
% Solve the problem
sol = optimize(Constraints,Objective,ops);
if sol.problem == 0
    % Extract and display value
	[row,col,v] = find(value(x)==1);
    MINLP_output = struct('i',row','j',col','m',meta_part,'c',value(y(meta_part)));
else
    display('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
end