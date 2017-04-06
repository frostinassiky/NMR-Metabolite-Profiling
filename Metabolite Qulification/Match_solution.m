%% Step 1 Data structure
% h_H = waitbar(0,'running H');
gauss_p1 = 0.2;gauss_p2 = 0.1;
x_local = old1_x;
clear single_link
lb1 = zeros(settings.num,1);ub1 = settings.conc_max*ones(settings.num,1);
lb2 = -0.03*ones(sum([metabolites.num]),1);
ub2 = -1*lb2;
%% A
[pks_s,locs_s] =findpeaks(spectrum,ppm,'MinPeakHeight' ,0.005*max(spectrum),'MinPeakProminence',0.005);%,'MinPeakWidth',0.001);
num_peak_s = length(locs_s);

%%  solve linear function array
A = zeros(num_peak_s);
b = zeros(num_peak_s,1);
for l = 1:num_peak_s
    b(l,1) = pks_s(l);
    for m = 1:num_peak_s
        A(l,m) = lorezf(locs_s(l),locs_s(m),1);
    end
end
pks_s = A\b;


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
[locs_th, index_locs] = sort(locs_th);
pks_th = pks_th(index_locs);
peak_index = peak_index(index_locs);


%% O
conc = zeros(settings.num,3);

%% C 
strength = zeros(settings.num_peak,num_peak_s);
%% Answer
answer = struct('name',[],'concentration',[],'lb',[],'ub',[]);

%% Step 2 Calculation
% conc
conc(:,1) = old1_x(1:settings.num);
conc(:,2) = ci(1:settings.num,1);
conc(:,3) = ci(1:settings.num,2);
% strength
for k = 1:settings.num_peak
for l = 1:num_peak_s
	if abs(locs_th(k)-locs_s(l))>0.03
		continue
	end
	strength(k,l) = gaussmf(locs_th(k),[gauss_p1, locs_s(l)]) * ...
		gaussmf(pks_th(k)*conc(peak_index(k)),[gauss_p2, pks_s(l)]);
end
end

%% Step 3 Involution (only Empty and Single links )
% waitbar(0.1, h_H, 'running H: Involution');
tot = sum(sum(strength>0));
involution = 1;
while involution
%     waitbar(0.1 + 0.9 - 0.9* sum(sum(strength>0))/tot , h_H, 'running H: Involution');
	% empty_link
    disp(sum(sum(strength>0)));
	empty_link = find( max(strength')==0 );
	for peak_del = empty_link
		% peak_del = empty_link(1);
		% change bound
        if pks_th(peak_del)>0.01&&...
                (~exist('single_link') || sum(max(strength(k,:)')==0) > length(k)/2 )
            lb1(peak_index(peak_del)) = 0;ub1(peak_index(peak_del)) = eps;
            conc(peak_index(peak_del)) = 0;
        	k = find( peak_index == peak_index(peak_del) );
   			strength(k,:) = -1 * ones( size(strength(k,:)) );
		end
		% empty_link = find( max(strength')==0 );
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
            lsqnonlin(fun_local,[conc(meta_cert,1),locs_s(peak_cert_s)-chemshift],...
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
				strength(k,l) = gaussmf(locs_th(k),[gauss_p1, locs_s(l)]) * ...
					gaussmf(pks_th(k)*conc(peak_index(k)),[gauss_p2, pks_s(l)]);
				% get best score
				if best_score < strength(k,l)
					best_score = strength(k,l);
					best_link = [k,l];
				end
			end
			% record and delete
            if length(best_link) == 2
                lb2(index_locs(k)) = locs_s(best_link(2)) - locs_th(best_link(1)) ;
                ub2(index_locs(k)) = lb2(index_locs(k))+eps;
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
    %% Partitioned Matrix
    % check xone_link less than 10
    xone = find(pks_s'>settings.pks_min*5 & max(strength)>0);
    x_min = 10;
    x_index = 0;
    for k = 1:length(xone)
        x = length( find(strength(:,xone(k))>0) );
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
    if xone_link
        [~,index] = max(pks_s(xone(x_index)));
        % the highest peaks in sle
		peak_xone_s = xone(x_index(index));
        % corresponding peaks and metabolites
        peak_xone_th = find(strength(:,peak_xone_s)>0);
		meta_xone = unique( peak_index( peak_xone_th ) );
        % relative peaks in theory
        peak_part_th = [];
        for k = meta_xone'
            peak_part_th = [peak_part_th;find( peak_index == k )];
        end
        % relative peaks in sle
        [~,peak_part_s] = find(strength(peak_part_th,:)>0);
        peak_part_s = unique(peak_part_s);
        % Partitioned Matrix
        pstrength = strength(peak_part_th,peak_part_s);
        %------------%
        eval( settings.minlp_solver )
        %------------%
        for k = length(MINLP_output.m)
            if find( peak_index(peak_part_th(MINLP_output.i)) == MINLP_output.m(k) )
                conc(MINLP_output.m(k),1) = MINLP_output.c(k);
            else
                conc(MINLP_output.m(k),1) = 0;
            end
        end
		% change bound
		lb1(MINLP_output.m) = conc(MINLP_output.m(k),1);
        ub1(MINLP_output.m) = conc(MINLP_output.m(k),1)+eps;
		% record and delete
        strength(peak_part_th,:) = -1;
        for k = 1:length(MINLP_output.i)
            best_link = ...
                [peak_part_th(MINLP_output.i(k)),peak_part_s(MINLP_output.j(k))];
            strength(:, best_link(2) ) = -1;
            strength( best_link(1), best_link(2) ) = -2;
            lb2(index_locs(best_link(1))) = locs_s(best_link(2)) - locs_th(best_link(1));
            ub2(index_locs(best_link(1))) = lb2(index_locs(best_link(1))) + eps;
        end
        surf(strength,'EdgeColor','None');view(2);axis equal;

		% re-find xone type link
%         UpdateConcAndStrength;
%         xone = find(pks_s'>settings.pks_min*5 & max(strength)>0);
%         x_min = 10;
%         x_index = 0;
%         for k = 1:length(xone)
%             x = length( find(strength(:,xone(k))>0) );
%             if x < x_min && x > 0
%                 x_min = x;
%                 x_index = k;
%             elseif x == x_min && x > 0
%                 x_index(end+1) = k;    
%             end
%         end
%         xone_link = 0;
%         if x_min < 10
%             xone_link = 1;
%         end
    end
    involution = sum(sum(strength>0));
end
%% Step 4 ready to unlinearfit
old2_x2 = zeros(settings. num_clust,1);
counter = 0;
counter_peak = 0;
for k = 1:length(metabolites)
    num_peak = metabolites(k).num;
    for l = 1:num_peak
        index = counter+metabolites(k).peak(l,3);
        factor = 1/length(find( metabolites(k).peak(l,3)==metabolites(k).peak(:,3)));
        old2_x2(index) = old2_x2(index) + ...
            (lb2(counter_peak+l)+ub2((counter_peak+l)))*factor/2;
    end
    
%         for l = 1:num_peak
%         index = counter+metabolites(k).peak(l,3);
% %         factor = 1/length(find( metabolites(k).peak(l,3)==metabolites(k).peak(:,3)));
%         factor = 0;
%         peak_bro = find( metabolites(k).peak(l,3)==metabolites(k).peak(:,3));
%         for m = 1:peak_bro
%             right_rat = 1/length(find( abs( metabolites(k).peak(m,2)- locs_s)<0.02));
%             factor = factor + right_rat;
%         end
%         right_rat = 1/length(find( abs( metabolites(k).peak(l,2)- locs_s)<0.02));
%         old2_x2(index) = old2_x2(index) + ...
%             (lb2(counter_peak+l)+ub2((counter_peak+l)))*right_rat/factor/2;
% 
%         if isnan(old2_x2(index) )
%             old2_x2(index) = 0;
%         end
%         end
    
        
    counter = counter + length(metabolites(k).clust);
    counter_peak = counter_peak + metabolites(k).num;
end

old2_x = [lb1;old2_x2];%-5*(ppm(2)-ppm(1))
x_save = old2_x(1:17);
time = fix(clock);
save(['old2_x_',num2str(time(5)),'-',num2str(time(6)),'.txt'],'x_save','-ascii');
