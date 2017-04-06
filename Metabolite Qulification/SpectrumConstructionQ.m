function [ppm, spectrum] = SpectrumConstructionQ(  concentration_minishift, para   )
%% SpectrumConstruction Construct Spectrum 
%   concentration_minishift is combined to make a better called
    metabolites = para.metabolites;
    num_metabolite = length(metabolites);
    conc = concentration_minishift(1:num_metabolite);
    index_minishift = ceil( concentration_minishift(num_metabolite+1:end) );
    minishift = zeros(size(index_minishift));
    for k = 1:length(index_minishift)
        try
            minishift(k) = para.set_minishift{k}(index_minishift(k));
        catch
            disp(k)
            disp(index_minishift(k))
            disp(para.set_minishift{k}(index_minishift(k)))
        end
    end
    ppm_sample1 = zeros(size(minishift));
    counter = 0;
    for k = 1:num_metabolite
        num_peak = metabolites(k).num;
        for l = 1:num_peak
            counter = counter + 1;
            ppm_sample1(counter) = metabolites(k).peak(l,2);
        end
    end
    ppm_sample2 = ppm_sample1 + minishift;
    ppm_sample = sort([ppm_sample1;ppm_sample2]);
    
    spectrum = zeros( size(ppm_sample) );
    counter = 0;
    for k = 1:num_metabolite
        spect = zeros( size(ppm_sample) );
        num_peak = metabolites(k).num;
        for l = 1:num_peak
            counter = counter + 1;
            lorez = lorezf(ppm_sample,...
                metabolites(k).peak(l,2)+minishift(counter),...
                metabolites(k).peak(l,1));
            % figure(1)
            % plot(ppm,lorez)
            spect = spect + lorez;
        end
        spectrum = spectrum + conc(k)*spect;  
    end
    ppm = ppm_sample;
end

