function [I,C] = index_set_I(I_delta)
    C=zeros(size(I_delta,1),1);
    idx_save=[];
    for i=1:1:size(I_delta,1)
        idx=I_delta(i,:);
        C(i)=c_coeff(I_delta,idx);
        if C(i)~=0
            idx_save=[idx_save,i];
        end
    end
    C=C(idx_save,:);
    I=I_delta(idx_save,:);
end