function next_earthquake_plot(states,repeats,settings)
    states=cellfun(@(state) state.time_after,states);
    [f,x,~,~] = ecdf(states,'frequency',repeats);
    hold off
    stairs(x(2:end),100*f(2:end),'linewidth',1.5);
    set(get(gca,'children'),'Color',[0.0,0.0,0.0])
    hold on
    set(gca,'XScale','log');
    xlim([0,settings.cutoff-eps]);
    grid on
    xlabel('Time until next event (yr)','interpreter','latex');
    ylabel('Probability (\%)','interpreter','latex')
end
