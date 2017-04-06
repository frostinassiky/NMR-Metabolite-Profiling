%% prepare YASoverf
%move the common task here
meta_part = unique(peak_index(peak_part_th));
        if sum(meta_part == 15)
            disp('here')
        end
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
% Set some options for YALMIP and solver
ops = sdpsettings('verbose',0,'solver','bnb','showprogress',0);

grid_fun = @(y) ...
    YASolverf_quick(x,y,Constraints,peak_part_th,peak_part_s,gauss_p1, gauss_p2, pks_th, locs_th,  pks_s,locs_s, ops, pstrength, peak_index,meta_order );

% close all;
try
	close(h)
end
%% Start parallel pool (parpool) using the 'local' profile
% while 1
%     p = gcp('nocreate'); % If no pool, do not create new one.
%     if isempty(p)
%         parpool
%     else
%         poolsize = p.NumWorkers;
%         break
%     end
% end   
%% Problem Definition

x0 = unifrnd(0,9*ones(size(meta_part)));
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','TolX',0.01);
bestx = fmincon(grid_fun,x0,[],[],[],[],0,9,[],opts);



%% output
[~,pairs] = grid_fun(bestx);
[row,col,v] = find(pairs>=0.95);
MINLP_output = struct('i',row','j',col','m',meta_part','c',bestx');



