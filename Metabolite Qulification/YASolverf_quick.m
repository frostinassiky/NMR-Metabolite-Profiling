function [score,pairs] = YASolverf_quick(x,y,Constraints,peak_part_th,peak_part_s,gauss_p1, gauss_p2, pks_th, locs_th,  pks_s,locs_s, ops, pstrength, peak_index,meta_order )
% change objective
    for k = 1:length(peak_part_th)
    for l = 1:length(peak_part_s)
        if pstrength(k,l) == 0
            continue;
        else
            pstrength(k,l) = gaussmf(locs_th(peak_part_th(k)),[gauss_p1, locs_s(peak_part_s(l))]) * ...
                gaussmf(pks_th(peak_part_th(k))*y(meta_order(peak_index(peak_part_th(k)))),[gauss_p2, pks_s(peak_part_s(l))]);
        end
    end
    end
    
    Objective = -sum(x(:).*pstrength(:));
   
    % Solve the problem
    optimize(Constraints,Objective,ops);
    score = value(Objective);
    pairs = value(x).*pstrength;
    
