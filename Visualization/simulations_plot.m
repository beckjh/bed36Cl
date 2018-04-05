function simulations_plot(states,fault,repeats,i_max,confidence)
    dg = [ 0.3 0.6 1];
    h0 = errorbar(fault.samples.h',fault.samples.y,...
        fault.samples.sigma','.','MarkerSize',20,...
        'MarkerEdgeColor',dg,'MarkerFaceColor',dg,'LineWidth',1.2);
    hold on
    A = cellfun(@(x) {x.simulation'},states,'UniformOutput',true);
    A = cell2mat(A);
    q=w_quantile(A,repeats,...
        [(confidence+(100-confidence)/2)/100,0.5,(100-confidence)/2/100]);
    h2 = plot(fault.samples.h,q(2,:),'r-');
    h2.Color = [0,0,0];
    h2.LineWidth = 2;
    hold on
    h = plot(fault.samples.h,q(1,:));
    h.Color = [1,0,0];
    h = plot(fault.samples.h,q(3,:));
    h.Color = [1,0,0];
    h3 = plot(fault.samples.h,A(i_max,:),'b--');
    h3.Color = [0,0,0];
    h3.LineWidth = 2;
    xlabel('Height (cm)','interpreter','latex')
    ylabel('Concentration of ${}^{36}$Cl (at/g)','interpreter','latex')
    legend([h0,h2,h,h3],{'Measurements','Medians',sprintf('$%.0f %s$-confidence bands',confidence,'\%'),'Least squares solution'},'Location','NorthWest','interpreter','latex');
end
