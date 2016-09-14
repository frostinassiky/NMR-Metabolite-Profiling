%% NMR spectrum simulation
% This is simultion of NMR spectrum, which based on Lorenz function
%
% By Frost @ UAlberta

%% Parameters
% the difault file is 
%%
% 
% * metabolites.xlsx
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
h = waitbar(0,'By Frost @ UAlberta');
waitbar(0.1,h,'Step1, Refresh database')
try 
    ! python database_update.py
catch 
    disp('Error, please add python in system path')
end

%% Step2, read settings
% There should be a setting file named settings.txt
waitbar(0.5,h,'Step1, Refresh database')
fileID = fopen('settings.txt');
settings = struct;
file = textscan(fileID,'%s','delimiter','\n'); 
for k = 1:length([file{1}])-2
    file_line{k} = file{1}{k+2};
end

global settings;
temp = textscan(file_line{1},'%s','delimiter','='); 
ppm_range = textscan(temp{1}{2},'%s','delimiter',','); 
settings.ppm_range = [str2double(ppm_range{1}{1}),str2double(ppm_range{1}{2})];

temp = textscan(file_line{2},'%s','delimiter','='); 
ppm_prec = textscan(temp{1}{2},'%s','delimiter',','); 
settings.ppm_prec = [str2double(ppm_prec{1}{1}),str2double(ppm_prec{1}{2})];

temp = textscan(file_line{3},'%s','delimiter','='); 
ppm_prec = textscan(temp{1}{2},'%s','delimiter',','); 
settings.var = str2double(ppm_prec{1}{1});
settings.pks_min = 0.001;

fclose(fileID);
%% 
% Then read _metabolites.xlsx_ to get concentration and shift
%
% They will all be saved in *settings*
[~,~,metabolites] = xlsread('metabolites.xlsx');
settings.name=metabolites(2:end,1);
settings.conc=metabolites(2:end,2);
settings.num = length(settings.name);
for k = 1:settings.num
    try
        shift=textscan(metabolites{k+1,3},'%s','delimiter',',');
        settings.shifts{k,1}=str2double(shift{1});
    catch
        settings.shifts{k,1}=[];
    end
end


%% Step3, spectum initial
% the length and precition depend on *settings* file
waitbar(0.75,h,'Step3, spectum initial')
if settings.ppm_prec(1)==1
    ppm = settings.ppm_range(1):settings.ppm_prec(2):settings.ppm_range(2);
else
    ppm = linspace(settings.ppm_range(1),settings.ppm_range(2),settings.ppm_prec(2));
end
spectrum = zeros(size(ppm));

%% Step4, merge each spectum 
% the function is 
% 
% $$S(ppm)=\Sigma \Sigma lorez(ppm-ppm_peak;concentrarion)$$
%%
% 
%   function res = lorezf(ppm,peak_ppm,peak_Height)
%       w=0.01;    
%       res = peak_Height*w./(w+4*(ppm-peak_ppm).^2);
%   end
% 
 

waitbar(0.80,h,'Step4, merge each spectum')
figure(1);hold on;set(figure(1),'visible','off')
figure(2);hold on;set(figure(2),'visible','off')
figure(3);hold on;set(figure(3),'visible','off')
list_random = [];
for k = 1:settings.num
    waitbar(0.80+0.2*k/settings.num,h,'Step4, merge each spectum..')
%     if settings.conc{k}==0
%         continue
%     end
	spect = zeros(size(ppm));

    metabolite  =  MetaPreparation( settings.name{k});
    peak_number = metabolite.num;
    clust_shift = zeros(size(metabolite.clust));
    % lorenz function    
    for l = 1:length(metabolite.clust)
        a = randn*settings.var;
        while(abs(a)>0.03)
            a = randn*settings.var;
        end
        a = roundn(a,-4);
        clust_shift(l) = a;
        list_random(end+1) = a;
    end
    for l = 1:peak_number
        shift_local = metabolite.peak(l,2)+clust_shift(metabolite.peak(l,3));
        lorez = lorezf(ppm,shift_local,metabolite.peak(l,1));
        figure(1)
        plot(ppm,lorez)
        spect = spect + lorez;
    end
    figure(2)
    plot(ppm,spect)
    spectrum = spectrum + settings.conc{k}*spect;
        
end

%% Step5, figure
% * figure 1 shows us every peaks
% * figure 2 shows us each single spectrum
% * figure 3 shows us NMR spectrum
figure(1);grid on;hold off
xlabel('ppm');set(gca,'xdir','reverse')
ylabel('absolute intensity')
title('Peaks (simulition)')
figure(2);grid on
xlabel('ppm');set(gca,'xdir','reverse')
ylabel('absolute intensity')
title('Single sprectrum (simulition)')
%legend(settings.name)
plot(ppm,-spectrum);hold off
figure(3)
plot(ppm,spectrum+0.01*randn(size(ppm)),'LineWidth',1.5)
plot(ppm,spectrum)
grid on;hold off
xlabel('ppm');set(gca,'xdir','reverse')
ylabel('absolute intensity')
title('NMR sprectrum (simulition)')
legend('with noise','precise')
%% 
% display three figures

set(figure(1),'visible','on')
set(figure(2),'visible','on')
set(figure(3),'visible','on')
A=[ppm',spectrum'+0.002*randn(size(ppm'))];
list_random = list_random';
save NMRSpectrum.txt A '-ascii'
save randlist.txt list_random '-ascii'
close(h)
savefig(figure(1),'figure(1)')
savefig(figure(2),'figure(2)')
savefig(figure(3),'figure(3)')
close(figure(3))
close(figure(2))
close(figure(1))
	