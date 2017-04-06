set_minishift = cell(settings.num_peak,1);
size_minishift = zeros(settings.num_peak,1);
p = [1,1];
meta_delete = [];
meta_certain = [];
for k = 1:settings.num_peak
    chemshift = metabolites( p(1) ).peak( p(2),2 );
    chemshift = roundn(chemshift,-3);
    index =  find(EQ(ppm,chemshift-0.03)) : find(EQ(ppm,chemshift+0.03));
    [pks,locs] = ...
        findpeaks(spectrum(index),ppm(index),'MinPeakHeight' ,0.02*max(spectrum(index)),'MinPeakProminence',0.005);
	if isempty(pks) && (isempty(meta_delete) || meta_delete(end)~=p(1))
		meta_delete(end+1) = p(1);
	end
    if length(pks)==1 && (isempty(meta_certain) || meta_certain(end,1)~=p(1))
        meta_certain(end+1,:) = [p(1),p(2)];
    end
		
    set_minishift{k} = locs-metabolites( p(1) ).peak( p(2),2 );
    size_minishift(k) = size(set_minishift{k},1);
    p(2) = p(2) + 1;
    try
        metabolites( p(1) ).peak( p(2),2 );
    catch
        p(2) = 1;
        p(1) = p(1) + 1;
    end
end
meta_delete = sort(meta_delete,'descend');
[meta_certain(:,1),index] = sort(meta_certain(:,1),'descend');
meta_certain(:,2) = meta_certain(index,2);
k1 = 1; k2 = 1;
% while 1
%     if k1 > length(meta_delete) & k2 > length(meta_certain)
%         break
%     end
%     
%     if k2 > length(meta_certain) | k1 <= length(meta_delete) & meta_delete(k1) >= meta_certain(k2)
%         settings.num = settings.num - 1;
%         settings.num_peak = settings.num_peak - metabolites( meta_delete(k1) ).num;
%         answer(end+1)=...
%             struct('name',settings.name{meta_delete(k1)},'concentration',0,'lb',0,'ub',0);
%         metabolites( meta_delete(k1) ) = [];
%         settings.name{ meta_delete(k1)} = [];
%         if k1<=length(meta_delete) && k2<=length(meta_certain) && meta_delete(k1) == meta_certain(k2,1)
%             k2 = k2 + 1;
%         end
%         k1 = k1 + 1;
%     end
%     if k1 > length(meta_delete) | k2 <= length(meta_certain) & meta_delete(k1) < meta_certain(k2)
%         
%         chemshift = metabolites( meta_certain(k2,1) ).peak( meta_certain(k2,2),2 );
%         chemshift = roundn(chemshift,-3);
%         conc_local = metabolites( meta_certain(k2,1) ).peak( meta_certain(k2,2),1 );
%         index =  find(EQ(ppm,chemshift-0.03)) : find(EQ(ppm,chemshift+0.03));
%         
%         options = optimoptions(@lsqnonlin,'Display','iter-detailed',...
%             'MaxFunEvals', settings.ites(1),'TolFun',1e-14);
%         fun_local = @(x) ...
%             spectrum(index) - lorezf(ppm(index),chemshift+x(2),x(1) * conc_local);
%         [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
%             lsqnonlin(fun_local,[settings.conc_max/2,0],[0,-0.03],[settings.conc_max,0.03],options);
%         ci = nlparci(x,residual,'jacobian',jacobian);
%         if (ci(1,2)-ci(1,1))/x(1) < 0.20
%             answer(end+1)=...
%                 struct('name',settings.name{meta_certain(k2,1)},'concentration',x(1),'lb',ci(1,1),'ub',ci(1,2));
%             settings.name{meta_certain(k2,1)} = [];
%             settings.num_peak = settings.num_peak - metabolites( meta_certain(k2,1) ).num;
%             settings.num = settings.num - 1;
%             metabolites( meta_certain(k2,1) ) = [];
%         elseif (ci(1,2)-ci(1,1))/x(1) > 10
%             settings.num = settings.num - 1;
%             settings.num_peak = settings.num_peak - metabolites( meta_certain(k2,1) ).num;
%             answer(end+1)=...
%                 struct('name',settings.name{meta_certain(k2,1)},'concentration',0,'lb',0,'ub',0);
%             metabolites( meta_certain(k2,1) ) = [];
%             settings.name{ meta_certain(k2,1)} = [];
%         end
%         k2 = k2 + 1;
%     end
% end
settings.name = settings.name(~cellfun(@isempty, settings.name));
set_minishift = cell(settings.num_peak,1);
size_minishift = zeros(settings.num_peak,1);
p = [1,1];
for k = 1:settings.num_peak
    chemshift = metabolites( p(1) ).peak( p(2),2 );
    chemshift = roundn(chemshift,-3);
    index =  find(EQ(ppm,chemshift-0.03)) : find(EQ(ppm,chemshift+0.03));
    [pks,locs] = ...
        findpeaks(spectrum(index),ppm(index),'MinPeakHeight' ,0.02*max(spectrum(index)),'MinPeakProminence',0.005);
    set_minishift{k} = locs-metabolites( p(1) ).peak( p(2),2 );
    size_minishift(k) = size(set_minishift{k},1);
    p(2) = p(2) + 1;
    try
        metabolites( p(1) ).peak( p(2),2 );
    catch
        p(2) = 1;
        p(1) = p(1) + 1;
    end
end


para.set_minishift=set_minishift;
para.size_minishift=size_minishift;
