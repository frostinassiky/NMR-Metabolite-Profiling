function [metabolite]  =  MetaPreparation( name_metabolite)
%% SpectumPreparation: prepare each spectum 

%% the function is 
% 
% $$S(ppm)=\Sigma \Sigma lorez(ppm-ppm_peak;concentrarion)$$
%%
% 
%   function res = lorezf(ppm,peak_ppm,peak_Height)
%       w=0.01;    
%       res = peak_Height*w./(w+4*(ppm-peak_ppm).^2);
%   end
% 
    if nargin == 0 
        name_metabolite = 'Alanine';
    end
    global settings
    
    fileID = fopen(['database\', name_metabolite, '.txt']);
    file = textscan(fileID,'%s','delimiter','\n'); 
    fclose(fileID);

    	%%  search Clusters
	%   search 'List of Multiplets'
    for clust_index=1:length(file{1})
        if(strcmp(file{1, 1}{clust_index, 1},'Table of Multiplets') || ...
           strcmp(file{1, 1}{clust_index, 1},'List of Multiplets')     )
            break
        end
    end
%  search first Multiplets
    while 1
        clust_index = clust_index + 1;
        check_clust = file{1}{clust_index, 1};
        if( isempty(check_clust) ) 
            continue    
        end
        if( check_clust(1) == '1' )  
            break   
        end
    end
%  save each Multiplets   
    clust_number = 0;
    clust_number_save = 0;
    while 1
        clust_info = textscan(file{1, 1}{clust_index+clust_number, 1},'%s','delimiter','\t'); 
        shift_local = str2double(clust_info{1}{2});
        clust_number = clust_number + 1;
        clust_number_save = clust_number_save +1;
        clust(clust_number_save) = shift_local;
        
        try
            if(strcmp(file{1, 1}{clust_index+clust_number, 1},''))
                break
            end 
        catch
            break
        end
    end

    %%  search Peaks
%   search 'Table of Peaks'
    for peaks_index=1:length(file{1})
        if(strcmp(file{1, 1}{peaks_index, 1},'Table of Peaks') || ...
           strcmp(file{1, 1}{peaks_index, 1},'List of Peaks')     )
            break
        end
    end
%  search first peak
    while 1
        peaks_index = peaks_index + 1;
        check_peak = file{1, 1}{peaks_index, 1};
        if( isempty(check_peak) ) 
            continue    
        end
        if( check_peak(1) == '1' )  
            break   
        end
    end
%  save each peaks    
    peak_number = 0;
    peak_number_save = 0;
    while 1
        peak_info = textscan(file{1, 1}{peaks_index+peak_number, 1},'%s','delimiter','\t'); 
        switch length(peak_info{1})
            case 3
                index_height = 3;
            case 4 
                index_height = 4;
            otherwise
                index_height = 6;
        end
        shift_local = str2double(peak_info{1}{2});
        hight_local = str2double(peak_info{1}{index_height});
        [~,clust_index_local] = min( abs(shift_local-clust) );
        if  hight_local < settings.pks_min || exist('peak','var') && shift_local == peak(peak_number_save).pos
            peak_number = peak_number + 1;
            if(strcmp(file{1, 1}{peaks_index+peak_number, 1},''))
                break
            end   
            continue
        end
        peak_number = peak_number + 1;
        peak_number_save = peak_number_save +1;
        peak(peak_number_save) = struct('pos',shift_local,'h',str2double(peak_info{1}{index_height}),'cluster',clust_index_local);
        
        if(strcmp(file{1, 1}{peaks_index+peak_number, 1},''))
            break
        end        
    end
	peak_number = peak_number_save;

	
%%  solve linear function array
    A = zeros(peak_number);
    b = zeros(peak_number,1);
    for l = 1:peak_number
        b(l,1) = peak(l).h;
        for m = 1:peak_number
            A(l,m) = lorezf(peak(l).pos,peak(m).pos,1);
        end
    end
    X = A\b;
%%  solve first order function array
    h_sam = zeros(peak_number,1);
    x_sam = zeros(peak_number,1);
    for l = 1:peak_number
        h_sam(l) = peak(l).h;
        x_sam(l) = peak(l).pos;
    end
    x0 = [X;x_sam];
    options = optimoptions('fsolve','Display','off','TolFun',1e-12);
    fun = @(x) sec_order_fit(x,h_sam,x_sam);
    if isempty(x0)
        disp(['no peaks in ',name_metabolite] )
        x = x0;
    else
        x = fsolve(fun,x0,options);
    end
%%  return value   
    num_peak = peak_number;
    metabolite = reshape(x,num_peak,2);
    % metabolite_concentration = metabolite(:,1);
    % metabolite_shift = metabolite(:,2);
    metabolite = struct('num',num_peak,'peak',[metabolite,[peak.cluster]'],'clust',clust);
        
end

