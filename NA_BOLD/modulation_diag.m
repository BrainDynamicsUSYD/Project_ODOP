function [visibilityd,MaxValued] = modulation_diag(matrix_didecr) 

% Modulation of profile
 
jja = 1;
for kka = 2:(length(matrix_didecr)-1)

    % find max's and store in new array
    store_pre_maxd = matrix_didecr(kka-1);
    store_pos_maxd = matrix_didecr(kka+1);
    if matrix_didecr(kka) > store_pre_maxd && matrix_didecr(kka) > store_pos_maxd
        store_maxd = matrix_didecr(kka);
        store_max1d(jja) = store_maxd;
        jja=jja+1;
    end
end

    % find min's and store in new array
store_min1d = zeros; 
jji = 1;
for kki=2:(length(matrix_didecr)-1)
    store_pre_mind = matrix_didecr(kki-1);
    store_pos_mind = matrix_didecr(kki+1);
    if matrix_didecr(kki) < store_pre_mind && matrix_didecr(kki) < store_pos_mind
        store_mind = matrix_didecr(kki);
        store_min1d(jji) = store_mind;
        jji=jji+1;
    end
end

% Max value and indice of each row of stored max
store_max1dt = store_max1d';
[MaxValued,indMaxValued] = max(store_max1dt);
% Getting the max neighour values (at the right) close to the max
kkm = 1;
for kkm=1:size(indMaxValued,2)
    if length(store_max1d) > indMaxValued 
        store_max1_posd = store_max1d(indMaxValued(kkm)+1);
        store_max1_pos1d(1,size(indMaxValued,2)) = store_max1_posd;
    elseif length(store_max1d) > 1
        store_max1_pos1d = store_max1d(indMaxValued(kkm)-1);
    else
        store_max1_pos1d = 0;
    end
end
% Max value and indice of each row of stored min
store_min1dt = store_min1d';
[~,indMinValued] = max(store_min1dt);
% Getting the minimum value close to the max
kkmi = 1;
for kkmi=1:size(store_min1d,1)
    if length(store_min1d) >= indMinValued
        store_min1_mind = store_min1d(indMinValued(kkmi));
        store_min1_min1d(1,size(indMaxValued,2)) = store_min1_mind;
    else
        store_min1_min1d = 0;
    end
end

[Maximumd,ind_maxd]=max(store_max1d(:));
% modulation related to maximum value
jjma = 1;
for kkma=2:length(store_max1d)
    if kkma <= ind_maxd && store_max1d(kkma) ~= 0 && store_max1d(kkma-1) ~= 0
       modul_maxd = 100-((store_max1d(kkma-1)*100)/store_max1d(kkma));
       modul1_maxd(jjma) = modul_maxd;
       jjma=jjma+1;        
    end
end

% Visibility:
% ((0.5*(peak1+peak2))-valey1)/((0.5*(peak1+peak2))+valey1);
visibilityd = zeros(1,size(indMaxValued,2));
for vv=1:size(indMaxValued,2)
    visibilityd(vv) = ((0.5*(MaxValued(vv)+store_max1_pos1d(vv)))-store_min1_min1d(vv))/((0.5*(MaxValued(vv)+store_max1_pos1d(vv)))+store_min1_min1d(vv));
end