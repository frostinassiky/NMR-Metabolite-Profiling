%% Step 1 Data structure
% h_H = waitbar(0,'running H');
gauss_p1 = 0.05;gauss_p2 = 0.1;
x_local = old1_x;
%% A
[pks_samp,locs_samp] =findpeaks(spectrum,ppm,'MinPeakHeight' ,0.005*max(spectrum),'MinPeakProminence',0.005,'MinPeakWidth',0.001);
num_peak_samp = length(locs_samp);

%%  solve linear function array
A = zeros(num_peak_samp);
b = zeros(num_peak_samp,1);
for l = 1:num_peak_samp
    b(l,1) = pks_samp(l);
    for m = 1:num_peak_samp
        A(l,m) = lorezf(locs_samp(l),locs_samp(m),1);
    end
end
pks_samp = A\b;


%% B
% initial
pks_th = zeros(settings.num_peak,1);
locs_th = zeros(settings.num_peak,1);
peak_index = zeros(settings.num_peak,1);
p = [1,1];
% load
for k = 1:settings.num_peak
    pks_th(k) = metabolites( p(1) ).peak( p(2),1 );
	locs_th(k) = metabolites( p(1) ).peak( p(2),2 );
	peak_index(k) = p(1);
     p(2) = p(2) + 1;
    try
        metabolites( p(1) ).peak( p(2),2 );
    catch
        p(2) = 1;
        p(1) = p(1) + 1;
    end
end
% sort
[locs_th, index] = sort(locs_th);
pks_th = pks_th(index);
peak_index = peak_index(index);
% % discard peaks less than pks_min
% index = find(pks_th>stpks_min);
% locs_th = locs_th(index);
% pks_th = pks_th(index);
% peak_index = peak_index(index);
%% O
conc = zeros(settings.num,3);

%% C 
strength = zeros(settings.num_peak,num_peak_samp);
%% Answer
answer = struct('name',[],'concentration',[],'lb',[],'ub',[]);

%% Step 2 Calculation
% conc
conc(:,1) = old1_x(1:settings.num);
conc(:,2) = ci(1:settings.num,1);
conc(:,3) = ci(1:settings.num,2);
% strength
for k = 1:settings.num_peak
for l = 1:num_peak_samp
	if abs(locs_th(k)-locs_samp(l))>0.03
		continue
	end
	strength(k,l) = gaussmf(locs_th(k),[gauss_p1, locs_samp(l)]) * ...
		gaussmf(pks_th(k)*conc(peak_index(k)),[gauss_p2, pks_samp(l)]);
end
end

%% Step 3 Involution
% waitbar(0.1, h_H, 'running H: Involution');
tot = sum(sum(strength>0));
involution = 1;
while involution
%     waitbar(0.1 + 0.9 - 0.9* sum(sum(strength>0))/tot , h_H, 'running H: Involution');
	% empty_link
	empty_link = find( sum(strength==0,2)==size(strength,2) );
	while empty_link
		peak_del = empty_link(1);
		% change bound
		lb1(peak_index(peak_del)) = 0;ub1(peak_index(peak_del)) = eps;
		for k = find( peak_index == peak_index(peak_del) )
            
			strength(k,:) = -1 * ones( size(strength(k,:)) );
		end
		empty_link = find( sum(strength==0,2)==size(strength,2) );
	end 

	% single_link 
	single_link = find( sum(strength>0,1) == 1 );
	while single_link
		peak_cert_s = single_link(1);
        peak_cert_th = find( strength(:,peak_cert_s)>0 );
		if length(find( strength(peak_cert_th,:)>0 )) >1
            single_link(1) = [];
            continue
        end
		meta_cert = peak_index( peak_cert_th );
		% calculate concentration
		chemshift = locs_th(peak_cert_th);
		index =  find(EQ(ppm,chemshift-0.03)) : find(EQ(ppm,chemshift+0.03));
        fun_local = @(x) ...
            spectrum(index) - lorezf(ppm(index),chemshift+x(2),x(1) * pks_th(peak_cert_th));
        [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqnonlin(fun_local,[conc(meta_cert,1),locs_samp(peak_cert_s)-chemshift],...
            [0,-0.03],[settings.conc_max,0.03],options);
        ci2 = nlparci(x,residual,'jacobian',jacobian);
        % record result
%         disp('+++++++++++')
%         disp(meta_cert)
%         disp('-----------')
		conc(meta_cert,:) = [ x(1), ci2(1,1), ci2(1,2) ];
		% change bound
		lb1(meta_cert) = x(1);ub1(meta_cert) = x(1)+eps;
		% record best link and delete other relevant peaks
		for k = find( peak_index == peak_index(peak_cert_th) )' 
			best_link = 0;
			best_score = 0;
			for l = 1 : size(strength,2)
				% re-compute strength
				if strength(k,l) <= 0
					continue
                end
				strength(k,l) = gaussmf(locs_th(k),[gauss_p1, locs_samp(l)]) * ...
					gaussmf(pks_th(k)*conc(peak_index(k)),[gauss_p2, pks_samp(l)]);
				% get best score
				if best_score < strength(k,l)
					best_score = strength(k,l);
					best_link = [k,l];
				end
			end
			% record and delete
            if length(best_link) == 2
                lb2(k) = locs_samp(best_link(2)) - locs_th(best_link(1));
                ub2(k) = lb2(k)+eps;
                strength(best_link(1),:) = -1 * ones(1,size(strength,2));
                strength(:,best_link(2)) = -1 * ones(size(strength,1),1);
                strength(best_link(1),best_link(2)) = -2;
            else
                strength(k,:) = -0.5 * ones(1,size(strength,2));
            end
            surf(strength,'EdgeColor','None');view(2);axis equal;
            1;
		end
		single_link = find( sum(strength>0,1) == 1 );
    end
    %% 1-X X-1
    UpdateConcAndStrength;
    % check xone_link less than 10
    xone = find(pks_samp'>settings.pks_min*5 & max(strength)>0);
    x_min = 10;
    x_index = 0;
    for k = 1:length(xone)
        x = length( find(strength(:,xone(k))>settings.strth_min) );
        if x < x_min && x > 0
            x_min = x;
            x_index = k;
        elseif x == x_min && x > 0
            x_index(end+1) = k;    
        end
    end
    xone_link = 0;
    if x_min < 10
        xone_link = 1;
    end
    while xone_link
        [row,col,v] = find( strength(:,xone(x_index)) );
        [~,index] = max(v);
		peak_xone_s = xone(x_index(col(index)));
        peak_xone_th = row(index);
		meta_xone = peak_index( peak_xone_th );
		% calculate concentration --copy--
		% calculate concentration
		chemshift = locs_th(peak_xone_th);
		index =  find(EQ(ppm,chemshift-0.03)) : find(EQ(ppm,chemshift+0.03));
        fun_local = @(x) ...
            spectrum(index) - lorezf(ppm(index),chemshift+x(2),x(1) * pks_th(peak_xone_th));
        [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqnonlin(fun_local,[conc(meta_xone,1),locs_samp(peak_xone_s)-chemshift],...
            [0,-0.03],[settings.conc_max,0.03],options);
        ci2 = nlparci(x,residual,'jacobian',jacobian);
        % record result
        disp('+++++++++++')
        disp(meta_xone)
        disp('-----------')
		conc(meta_xone,:) = [ x(1), ci2(1,1), ci2(1,2) ];
		% change bound
		lb1(meta_xone) = x(1);ub1(meta_xone) = x(1)+eps;
		% record best link and delete other relevant peaks
		for k = find( peak_index == peak_index(peak_xone_th) )' 
			best_link = 0;
			best_score = 0;
			for l = 1 : size(strength,2)
				% re-compute strength
				if strength(k,l) <= 0
					continue
                end
				strength(k,l) = gaussmf(locs_th(k),[gauss_p1, locs_samp(l)]) * ...
					gaussmf(pks_th(k)*conc(peak_index(k)),[gauss_p2, pks_samp(l)]);
				% get best score
				if best_score < strength(k,l)
					best_score = strength(k,l);
					best_link = [k,l];
				end
			end
			% record and delete
            if length(best_link) == 2
                lb2(k) = locs_samp(best_link(2)) - locs_th(best_link(1));
                ub2(k) = lb2(k)+eps;
                strength(best_link(1),:) = -1 * ones(1,size(strength,2));
                strength(:,best_link(2)) = -1 * ones(size(strength,1),1);
                strength(best_link(1),best_link(2)) = -2;
            else
                strength(k,:) = -0.5 * ones(1,size(strength,2));
            end
            surf(strength,'EdgeColor','None');view(2);axis equal;
            1;
		end
		% re-find xone type link
        UpdateConcAndStrength;
        xone = find(pks_samp'>settings.pks_min*5 & max(strength)>0);
        x_min = 10;
        x_index = 0;
        for k = 1:length(xone)
            x = length( find(strength(:,xone(k))>settings.strth_min) );
            if x < x_min && x > 0
                x_min = x;
                x_index = k;
            elseif x == x_min && x > 0
                x_index(end+1) = k;    
            end
        end
        xone_link = 0;
        if x_min < 0;
            xone_link = 1;
        end
    end
%%        
        
        
        
        
        
	% if there is NO empty link or single link
	% just choose the best one
    num_choose = 2;
	if isempty(empty_link) && isempty(single_link) 
		% but firstly strength matrix should be updated
		UpdateConcAndStrength;
        
		[ k , l ] = find( strength==max(max(strength)) );
		if length(k)>10
			break
		end
		peak_best_th = k(1);
		peak_best_s = l(1);
		meta_best = peak_index( peak_best_th );
        
% 	while single_link
% 		peak_cert_s = single_link(1);
%         peak_cert_th = find( strength(:,peak_cert_s)>0 );
% 		if length(find( strength(peak_cert_th,:)>0 )) >1
%             single_link(1) = [];
%             continue
%         end
% 		meta_cert = peak_index( peak_cert_th );
        
        
        
        
		% then i copy the code from single link
		% and replace cert to best
		% calculate concentration
		chemshift = locs_th(peak_best_th);
		index =  find(EQ2(ppm,chemshift-0.03)) : find(EQ2(ppm,chemshift+0.03));
        options.TolFun=1e-14;
        fun_local = @(x) ...
            spectrum(index) - lorezf(ppm(index),chemshift+x(2),x(1) * pks_th(peak_best_th));
        [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqnonlin(fun_local,[conc(meta_best,1),locs_samp(peak_best_s)-chemshift],...
            [0,-0.03],[settings.conc_max,0.03],options);
        ci2 = nlparci(x,residual,'jacobian',jacobian);
        % record result
		conc(meta_best,:) = [ x(1), ci2(1,1), ci2(1,2) ];
		% change bound
		lb1(meta_best) = x(1);ub1(meta_best) = x(1)+eps;
		% record best link and delete other relevant peaks
		for k = find( peak_index == peak_index(peak_best_th) )' 
			best_link = 0;
			best_score = 0;
			for l = 1 : size(strength,2)
				% re-compute strength
				if strength(k,l) <= 0
					continue
				end
				strength(k,l) = gaussmf(locs_th(k),[gauss_p1, locs_th(k)]) * ...
					gaussmf(pks_th(k)*conc(peak_index(k)),[gauss_p2, pks_th(k)]);
				% get best score
				if best_score < strength(k,l)
					best_score = strength(k,l);
					best_link = [k,l];
				end
			end
			% record and delete
            if length(best_link) == 2
                lb2(k) = locs_samp(best_link(2)) - locs_th(best_link(1))-eps;
                ub2(k) = lb2(k)+2*eps;
                strength(best_link(1),:) = -1 * ones(1,size(strength,2));
                strength(:,best_link(2)) = -1 * ones(size(strength,1),1);
                strength(best_link(1),best_link(2)) = -2;
            else
                strength(k,:) = -0.5 * ones(1,size(strength,2));
            end
            surf(strength,'EdgeColor','None');view(2);axis equal;
            1;
		end
	end
end

%% Step 4 ready to unlinearfit
old2_x = [lb1;lb2];
