function spect =  SpectumPreparation( ppm,name_metabolite)
%% SpectumPreparation: prepare each spectum 
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
    spect = zeros(length(ppm),1 );
    fileID = fopen(['database\', name_metabolite, '.txt']);
    file = textscan(fileID,'%s','delimiter','\n'); 
    fclose(fileID);
%%  search 'Table of Peaks'
    for peaks_index=1:length(file{1})
        if(strcmp(file{1, 1}{peaks_index, 1},'Table of Peaks'))
            break
        end
    end
%%  search first peak
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
%%  save each peaks    
    peak_number = 0;
    while 1
        peak_info = textscan(file{1, 1}{peaks_index+peak_number, 1},'%s',...
            'delimiter','\t'); 
        shift_local = str2double(peak_info{1}{2,1});
        %shift = shift_local + shifts(peak_number+1); minishift should
        %occur in the next step
        if peak_number > 0 && shift_local == peak(peak_number).pos
            continue
        end
        peak_number = peak_number + 1;
        peak(peak_number) = struct('pos',shift,'h',str2double(peak_info{1}{4,1}));

        if(strcmp(file{1, 1}{peaks_index+peak_number, 1},''))
            break
        end        
    end

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
%%  solve second order function array
    h_sam = zeros(peak_number,1);
    x_sam = zeros(peak_number,1);
    for l = 1:peak_number
        h_sam(l) = peak(l).h;
        x_sam(l) = peak(l).pos;
    end
    x0 = [X;x_sam];
    options = optimoptions('fsolve','Display','','TolFun',1e-12);
    fun = @(x) sec_order_fit(x,h_sam,x_sam);
    x = fsolve(fun,x0,options);
%%  lorenz function    
    for l = 1:peak_number
        lorez = lorezf(ppm,x(l+peak_number)+shif,x(l));
        % figure(1)
        % plot(ppm,lorez)
        spect = spect + lorez;
    end
    % figure(2)
    % plot(ppm,spect)
end

