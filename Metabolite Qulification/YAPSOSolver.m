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

close all;
try
	close(h)
end
%% Start parallel pool (parpool) using the 'local' profile

if isempty(gcp('nocreate'))
    parpool
end
%% Problem Definition

CostFunction=@(x) pso_function(x);        % Cost Function

meta_part = unique(peak_index(peak_part_th));
nVar=length(meta_part);            % Number of Decision Variables

VarSize=[nVar 1];                  % Size of Decision Variables Matrix

VarMin=zeros(nVar,1);              % Lower Bound of Variables
VarMax=9*ones(nVar,1);             % Upper Bound of Variables
%% PSO Parameters

MaxIt=settings.ites(2);      % Maximum Number of Iterations

nPop=10;        % Population Size (Swarm Size)

% PSO Parameters
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient
% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % Inertia Weight
% wdamp=1;        % Inertia Weight Damping Ratio
% c1=chi*phi1;    % Personal Learning Coefficient
% c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.02*(VarMax-VarMin);
VelMin=-VelMax;

%%%
pso
%%%

disp(BestSol)
