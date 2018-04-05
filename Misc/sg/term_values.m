function [CYIM]=term_values(mapObj,func,I,C,knots)
    CYIM=[];
    for i=1:1:size(I,1)
        idx=I(i,:);
        M_idx=i2m(idx);
        %[M1_idx,M2_idx,M3_idx]=ndgrid(1:M_idx(1),1:M_idx(2),1:M_idx(3));
        I_new=multiidx_gen(length(idx),@(i,d) i(d)-1,[M_idx-1,max(M_idx)],1);
        for j=1:size(I_new,1)    
            x_j=[];x_j_str=[];
            for d=1:1:length(I_new(j,:))
                x_j=[x_j,knots(idx(d),I_new(j,d),d)];
                x_j_str=[x_j_str,sprintf('%.14f',x_j(d)),',']; % can remove last ,
            end
            if mapObj.isKey(x_j_str)
                y_j=mapObj(x_j_str);
            else
                y_j=func(x_j);
                mapObj(x_j_str)=y_j;              
            end
            CYIM=[CYIM;C(i),y_j,idx,I_new(j,:)];
        end
    end
end
