% test

%% prepare TASoverf
%move the common task here
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
% Set some options for YALMIP and solver
ops = sdpsettings('verbose',0,'solver','bnb','showprogress',0);

pso_function = @(y) ...
    YASolverf_quick(x,y,Constraints,peak_part_th,peak_part_s,gauss_p1, gauss_p2, pks_th, locs_th,  pks_s,locs_s, ops, pstrength, peak_index,meta_order );

e = 0.02;
result = zeros(9/e,1);
k=1;
for y=0:e:9
    result(k) = -pso_function(y);
    k = k+1;
    %isp(y)
end
testresult = figure;
plot(0:e:9,result,'LineWidth',2);
xlabel('concentration');
ylabel('result');
grid on;
savefig(testresult,'testresult')
delete(gcp('nocreate'))