function cv=CV(state,plot_time)
    times=state.times((state.times>plot_time(1))&(state.times<plot_time(2)));
    m=(plot_time(2)-plot_time(1))/length(times);
    s=sqrt(sum((times(2:end)-times(1:end-1)).^2)/(length(times)-1));
    cv=s/m;
end