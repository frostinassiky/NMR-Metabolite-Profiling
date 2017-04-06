function [spectrum,spect] = SpectrumConstruction( ppm, concentration_minishift, metabolites   )
%% SpectrumConstruction Construct Spectrum 
%   concentration_minishift is combined to make a better called
    num_metabolite = length(metabolites);
    conc = concentration_minishift(1:num_metabolite);
    minishift = concentration_minishift(num_metabolite+1:end);
    spectrum = zeros( size(ppm) );
    spect = zeros( size(ppm,1),num_metabolite );
    counter = 0;
    for k = 1:num_metabolite
        % spect = zeros( size(ppm) );
        num_peak = metabolites(k).num;
        for l = 1:num_peak
            shift_loc = metabolites(k).peak(l,2)+...
                minishift(counter+metabolites(k).peak(l,3));
            lorez = lorezf(ppm,shift_loc,...
                metabolites(k).peak(l,1));
            % figure(1)
            % plot(ppm,lorez)
            spect(:,k) = spect(:,k) + lorez;
        end
        counter = counter + length(metabolites(k).clust);
        spect(:,k) = spect(:,k)*conc(k);
        spectrum = spectrum + spect(:,k);  
    end

end

