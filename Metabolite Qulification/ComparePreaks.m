function score = ComparePreaks( spect1, spect2, n )
%ComparePreaks Summary of this function goes here
%   Detailed explanation goes here
n = 50;%floor(n+1);
[pks1,locs1] = findpeaks(spect1,'NPeaks',n,'SortStr','descend');
[pks2,locs2] = findpeaks(spect2,'NPeaks',n,'SortStr','descend');

% n = min(length(pks1),length(pks2));
% [pks1,locs1] = findpeaks(spect1,'NPeaks',n);
% [pks2,locs2] = findpeaks(spect2,'NPeaks',n);
try
    score = ( sort(locs1)- sort(locs2)) / ( norm(locs2,2) ); 
catch
    score =100000000000;
end

end

