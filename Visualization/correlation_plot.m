function MCMC_correlation_plot(states,repeats,func,labels)
    A=cellfun(func,states,'UniformOutput',false);
    A=cell2mat(A')';
    A=repelem(A,repeats,1);
    [h,ax,~,H]=plotmatrix(A,A);
    for i=1:size(A,2)
        ax(i,1).YLabel.String=labels{i};      
    end
    for i=1:size(A,2)
        ax(size(A,2),i).XLabel.String=labels{i};
    end
    for i=1:size(h,1)
        for j=1:size(h,2)
            %h(i,j).Color=[100,100,100]/255.0;
            if i==j
                %h(i,i).FaceColor=[150,150,150]/255.0;
            end
        end
    end
    for i=1:size(H,2)
        H(i).FaceColor=[100,100,100]/255.0;
    end
    %ax(3,1).XLabel.String='Test7'; 
    %ax(3,2).XLabel.String='Test8'; 
    %ax(3,3).XLabel.String='Test9';
end