function LS_plot(scenario,states,prior,settings,i_max,plot_time)
    sampler=Sampler(scenario,prior,[],settings);
    sampler.N_bm=length(states{i_max}.bm);
    [~,X]=sampler.find_times(states{i_max});
    times=linspace(settings.T_min,settings.T_max,length(X));
    i_plot=max(1,max_such_that(times,@(x) x<plot_time));
    plot(times(i_plot:end)/1000.0,X(i_plot:end),'r');
    hold on
    scatter(states{i_max}.times(2:end)/1000.0,zeros(1,length(states{i_max}.times(2:end))),max(states{i_max}.jumps(2:end),1),'r','filled')
    hold on 
    plot(states{i_max}.freq_switches/1000.0,zeros(1,length(states{i_max}.freq_switches)),'k.','MarkerSize',20);
    hold on
    xlabel('Time (kyr)','interpreter','latex')
    axis([plot_time(1)/1000.0,plot_time(2)/1000.0,0,1])
    hold off
    legend({'Brownian motion','Earthquakes','Frequency changes'},'Location','NorthWest')
end