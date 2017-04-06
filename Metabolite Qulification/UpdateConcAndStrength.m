% nonlim
options.TolFun=0.3;
[x,~,residual,~,~,~,jacobian] = ...
    lsqnonlin(fun_sle,[conc(:,1);x_local(settings.num+1:end)],...
    [lb1;lb2],[ub1;ub2],options);
ci3 = nlparci(x,residual,'jacobian',jacobian);
x_local = x;

% conc
conc(:,1) = x_local(1:settings.num);
conc(:,2) = ci3(1:settings.num,1);
conc(:,3) = ci3(1:settings.num,2);
% strength
for k = 1:settings.num_peak
for l = 1:num_peak_s
	if strength(k,l) <= 0
		continue
	end
	strength(k,l) = gaussmf(locs_th(k),[gauss_p1, locs_s(l)]) * ...
		gaussmf(pks_th(k)*conc(peak_index(k)),[gauss_p2, pks_s(l)]);
end
end