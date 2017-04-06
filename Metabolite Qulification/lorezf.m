function res = lorezf(ppm,peak_ppm,peak_Height)
	% the prior value of w is 1 Hz/F, where F = 500Hz
	% So w = 1/500
    % w=0.000008;   
    w=0.002;
    % res = peak_Height*w./(w+4*(ppm-peak_ppm).^2);
    w = w/2;
    res = peak_Height* w^2 ./(w^2+(ppm-peak_ppm).^2);
    % ref: http://mathworld.wolfram.com/LorentzianFunction.html
end