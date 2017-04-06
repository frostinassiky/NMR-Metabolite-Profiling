%% Metabolite Quantification
% Quantify the Metabolites from NMR spectrum, which based on LS.
%
% By Frost @ UAlberta
warning off
addpath(genpath('C:\Users\Frost\Documents\IBM\ILOG\CPLEX_Studio_Community1263\cplex\examples\src\matlab'))
addpath(genpath('E:\Program Files\IBM\ILOG\CPLEX_Studio_Community1263\cplex\matlab\x64_win64'))
addpath(genpath('yalmip'))
%% Parameters
% the difault file is 
%%
% 
% * metabolites list.xlsx
% * settings.txt
% * lorezf.m
% * database_update.py
%

%% Step1, Refresh database
% Call function _databse_update_ , this file is write by *python*
%
% It will not take too much if we have got the database
%
% Files download from <http://www.hmdb.ca here>
clc;clear
h = waitbar(0,'By Frost @ UAlberta'); tic;
waitbar(0.1,h,'Step1, Refresh database')
try 
    ! python database_update.py
catch 
    disp('Error, please add python in system path')
end


%% Step2, read settings
% There should be a setting file named settings.txt
waitbar(0.2,h,'Step2, read settings')
fileID = fopen('settings.txt');
global settings
settings = struct;
file = textscan(fileID,'%s','delimiter','\n'); 
for k = 1:length([file{1}])-2
    file_line{k} = file{1}{k+2};
end

temp = textscan(file_line{5},'%s','delimiter','='); 
ites = textscan(temp{1}{2},'%s','delimiter',','); 
settings.ites = ...
    [str2double(ites{1}{1}),str2double(ites{1}{2}),...
    str2double(ites{1}{3}),str2double(ites{1}{4})];

temp = textscan(file_line{6},'%s','delimiter','='); 
ites = textscan(temp{1}{2},'%s','delimiter',','); 
settings.optimal_method = ites{1}{1};

temp = textscan(file_line{7},'%s','delimiter','='); 
ites = textscan(temp{1}{2},'%s','delimiter',','); 
settings.conc_max = str2double(ites{1}{1});

temp = textscan(file_line{8},'%s','delimiter','='); 
ites = textscan(temp{1}{2},'%s','delimiter',','); 
settings.pks_min = str2double(ites{1}{1});

temp = textscan(file_line{9},'%s','delimiter','='); 
ites = textscan(temp{1}{2},'%s','delimiter',','); 
settings.strth_min = str2double(ites{1}{1});

temp = textscan(file_line{10},'%s','delimiter','='); 
ites = textscan(temp{1}{2},'%s','delimiter',','); 
settings.minlp_solver = ites{1}{1};

% temp = textscan(file_line{2},'%s','delimiter','='); 
% ppm_prec = textscan(temp{1}{2},'%s','delimiter',','); 
% settings.ppm_prec = [str2double(ppm_prec{1}{1}),str2double(ppm_prec{1}{2})];
% 
% fclose(fileID);

%% Step3, read metabolites list
% read _metabolites list.xlsx_ to get names
%
% They will all be saved in *settings*
waitbar(0.3,h,'Step3, read metabolites list')
[~,~,metabolites] = xlsread('metabolites list.xlsx');
settings.name=metabolites(2:end,1);
settings.num = length(settings.name);
metabolites = struct('num',{},'peak',{},'clust',{});
settings.num_peak = 0;
settings.num_clust = 0;
for k = 1:settings.num
    metabolites(k) = MetaPreparation( settings.name{k} );
    settings.num_peak = settings.num_peak + metabolites(k).num;
    settings.num_clust = settings.num_clust + length(metabolites(k).clust);
end

%% Step4, NMR spectum initial
% Our gaol is saved in NMRspectrum.txt
%
% * ppm_sample is a subset from ppm
% * spectrum_sample is a subset from spectrum
% 
waitbar(0.4,h,'Step4, NMR spectum initial')
A = dlmread('NMRspectrum.txt');
ppm = A(:,1);
spectrum = A(:,2);
index = spectrum> 0.10*max(spectrum);
index = imdilate( index, ones( ceil(0.03/(ppm(2)-ppm(1))), 1) );
ppm_sample = ppm( index ) ;
spectrum_sample = spectrum( index );

%% Step5, Prepare initial vector of sulotion
% X0 consist of 2 part
% the first part is Cencentration: settings.num * 1
% the last part is minishift: 
%%
% 
% $$\Sigma_{1}^{settings.num}shifts per metabolite$$
% 
waitbar(0.5,h,'Step5, Prepare initial vector of sulotion')
x1 = ones( settings. num,1 );
%% 
% X2 is much more important than X1
% 
% That is because _X1_ is nonlinear.
% And the FWHM is 0.02, less than our max *minishift*, 0.03.
%
% So firstly, convolution is made to find a proper chemshift centre.

x2 = zeros( settings. num_clust,1 );
% p = [1,1];
% e = ppm(2) - ppm(1);
% tol = eps(2); % reset equal
% for k = 1:settings.num_peak
%     chemshift = metabolites( p(1) ).peak( p(2),2 );
%     chemshift = roundn(chemshift,-3);
%     index =  find(EQ(ppm,chemshift-0.03)) : find(EQ(ppm,chemshift+0.03));
%     period = lorezf( -0.03:e:0.03,0,1);
%     temp = conv( spectrum(index), period);
%     % find local maximal large value
%     temp = find( (temp>0.25*max(temp)) .* imregionalmax(temp) );
%     % find neareat one
%     pos_max = abs(temp - length(index));
%     pos_best = find(EQ(pos_max, min(pos_max)),1);
%     try
%         x2(k) = ( temp(pos_best)- length(index))*e;
%     catch
%         x2(k) = 0;
%     end
%     p(2) = p(2) + 1;
%     try
%         metabolites( p(1) ).peak( p(2),2 );
%     catch
%         p(2) = 1;
%         p(1) = p(1) + 1;
%     end
% end
x0 = [x1;x2];

%% Step6, Prepare two bonds
waitbar(0.6,h,'Step6, two bonds')

lb1 = zeros(size(x1));
ub1 = settings.conc_max*ones(size(x1));
lb2 = -0.03*ones(size(x2));
ub2 = 0.03*ones(size(x2));
lb = [lb1;lb2];
ub = [ub1;ub2];    

%% Step7, Solve nonlinear least-squares (nonlinear data-fitting) problems
% Both ppm and ppm_sample are availible
%

%% 
% First step, optimal concentration, only sample data is used
% 

waitbar(0.7,h,'Step7, Solve nonlinear least-squares problems')
options = optimoptions(@lsqnonlin,'Display','iter-detailed',...
    'MaxFunEvals', settings.ites(1),'TolFun',1e-14);
fun_sample = @(x) ...
    SpectrumConstruction( ppm_sample, x, metabolites )-spectrum_sample;
% fun = @(x) ...
%     SpectrumConstruction( ppm, x, metabolites )-spectrum;
% fun = @(x) ...
%     ComparePreaks( SpectrumConstruction( ppm, x, metabolites )...
%     , spectrum, length(x)/2 );

[x,~,residual,~,~,~,jacobian] = lsqnonlin(fun_sample,x0,lb,ub,options);
ci = nlparci(x,residual,'jacobian',jacobian);

old1_x = x;


%% 
% Then run settings.optimal_method
% old2_x is set in eval
waitbar(0.7,h,'Step7.3, run settings.optimal\_method')
spectrum_copy=spectrum;
spectrum=smooth(spectrum,5);
eval(settings.optimal_method);
spectrum=spectrum_copy;

% index_minishift = ceil( BestSol.Position );
% minishift = zeros(size(index_minishift));
% for k = 1:length(index_minishift)
%     minishift(k) = para.set_minishift{k}(index_minishift(k));
% end
%  
% old2_x = [old1_x(1:settings.num);minishift];

%% 
% At last, global optimation start
% 
h = waitbar(0.7, 'Step7.6, Last Step: Re-solve nonlinear least-squares problems');

% % Prepare two bonds
lb1 = zeros(size(x1));
ub1 = settings.conc_max*ones(size(x1));
for k = 1:settings.num
    if(old2_x(k)==0)
        ub1(k) = eps(1);
    end
end

lb2 = -0.03*ones(size(x2));
ub2 = 0.03*ones(size(x2));
lb = [lb1;lb2];
ub = [ub1;ub2];    

fun = @(x) ...
    SpectrumConstruction( ppm, x, metabolites )-spectrum;
old3_x = old2_x;
% for k = settings.num+1 :length(old2_x);
%     subfun = @(x) ...
%         norm(fun([old3_x(1:k-1)',x,old3_x(k:length(old2_x))']'),2);
%     old3_x(k) = fmincon(subfun,old2_x(k),[],[],[],[],-0.03,+0.03); 
% end    
%  


% for k = settings.num+1 :length(old2_x);
%     subfun = @(x) ...
%         norm(fun([old3_x(1:k-1)',x,old3_x(k:length(old2_x))']'),2);
%     test_x = -0.03;0.0001:0.03;
%     test_res = zeros(length(test_x),1);
%     for l = 1:length(test_x)
%         test_res(l) = subfun(test_x(l));
%     end
%     [~,index] = min(test_res);
%     old3_x(k) = fmincon(subfun,test_x(index),[],[],[],[],-0.03,+0.03); 
% end   

options = optimoptions(@lsqnonlin,'Display','iter-detailed',...
    'MaxFunEvals', settings.ites(3),'TolX',1e-4);
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(fun,old3_x,lb,ub,options); 

% 
% %% Step8, Nonlinear regression parameter confidence intervals
waitbar(0.8,h,'Step8, Nonlinear regression parameter confidence intervals')
% ci = nlparci(x,residual,'jacobian',jacobian);
% 
% %% Step9, Show result
waitbar(0.9,h,'Step9, Show result')
figure(2);grid on;hold on
xlabel('ppm');set(gca,'xdir','reverse')
ylabel('absolute intensity')
title(['Comparation (score) ']);%, num2str( norm(fun(x),2) )])

[cspectrum,spect] = SpectrumConstruction( ppm, x, metabolites );
plot(ppm,-spect,'LineWidth',0.5);

plot(ppm,spectrum,'LineWidth',1.5)
plot(ppm,cspectrum,'LineWidth',1)
legend(settings.name)%,'Sample','SpectrumConstruction')

% legend('Sample','SpectrumConstruction')
% disp('ci = nlparci(x,residual,_jacobian_,jacobian);')
% disp(ci)
close(h)
toc

