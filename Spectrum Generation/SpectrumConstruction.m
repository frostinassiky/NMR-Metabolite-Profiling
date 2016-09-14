function spectrum = SpectrumConstruction( ppm, concentration_minishift, metabolites   )
%SpectrumConstruction Construct Spectrum 
%   concentration_minishift is combined to make a better called
    num_metabolite = length(metabolites);
    conc = concentration_minishift(1:num_metabolite);
    minishift = concentration_minishift(num_metabolite+1:end);
    spectrum = zeros( size(ppm) );
    counter = 0;
    for k = 1:num_metabolite
        spect = zeros( size(ppm) );
        num_peak = metabolites(k).num;
        for l = 1:num_peak
            counter = counter + 1;
            lorez = lorezf(ppm,...
                metabolites(k).peak(l,2)+minishift(counter),...
                metabolites(k).peak(l,1));
            % figure(1)
            % plot(ppm,lorez)
            spect = spect + lorez;
        end
        spectrum = spectrum + conc(k)*spect;  
    end

end

