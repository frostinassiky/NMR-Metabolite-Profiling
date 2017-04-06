function score = ComparePreaks2( x, para, ppm, spectrum )
%ComparePreaks Summary of this function goes here
%   Detailed explanation goes here
    [ppm_sample, spect1] = SpectrumConstructionQ( x, para );
    index = zeros(size(spect1));
    for k = 1 : length(index)
        index(k) = find( EQ2(ppm,ppm_sample(k)) );
    end
    score = spect1 - spectrum(index);

end

