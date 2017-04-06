function y = fun( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    load('specrtum.mat')
    %y = SpectrumConstruction( ppm, x, metabolites )-spectrum;
    ComparePreaks( SpectrumConstruction( ppm, x, metabolites ), spectrum, n )
    
end

