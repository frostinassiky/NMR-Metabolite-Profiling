%% prepare YASoverf
%move the common task here
meta_part = unique(peak_index(peak_part_th));
        if sum(meta_part == 17)
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
switch length(meta_part)
    case 1                
		e = 5;
		result = zeros(settings.conc_max*e+1,1);
		x = 0:1/e:settings.conc_max;
		for k=1:length(x)
			result(k) = grid_fun(x(k));   
		end

		%% fnimcon opt
		opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','TolX',0.01*settings.conc_max);
		[~,index] = min(result);
		bestx = fmincon(grid_fun,x(index),[],[],[],[],0,settings.conc_max,[],opts);

    otherwise
        problem = createOptimProblem('fmincon',	'x0',5*ones(size(meta_part)),'objective',grid_fun,'options',opts);
		gs = GlobalSearch('Display','iter');
		bestx = run(gs,problem);

end
%% output
[~,pairs] = grid_fun(bestx);
[row,col,v] = find(pairs>=0.7);
MINLP_output = struct('i',row','j',col','m',meta_part','c',bestx');
