function [C] = c_coeff(index_set,idx)
    C=0;
    for i=1:1:size(index_set,1)
        idx_i=index_set(i,:);
        add=1;
        absi=0.0;
        for l=1:1:length(idx)
            if (idx_i(l)-idx(l))<0 || (idx_i(l)-idx(l))>1,
                add=0;
                break;
            end
            absi=absi+abs(idx_i(l)-idx(l));
        end 
        if add==1
            C=C+(-1)^absi;
        end
    end
end
