function BED_plots(results,p_burn_in,plot_time,confidence,range,resolution_factor)
    N_figs = 8; 
    n_figure = 1;
    clear_figs()
    pause(0.01)
    [states,repeats] = burn_in(results,p_burn_in,range);
    settings = results.settings;
    fault = results.fault;
    states = states(repeats>0);
    repeats = repeats(repeats>0);
    fprintf('#######################################################################################\n');
    fprintf('   %d MCMC iterations ## Burn in: %.1f%s ## Time window: %dka--%dka  \n',sum(results.n_proposed{1}),100*p_burn_in,'%',-plot_time(1)/1000,-plot_time(2)/1000);
    fprintf('#######################################################################################\n');     
    if sum(repeats)>0
        wRMSs = wRMS(cell2mat(cellfun(@(state) state.simulation',states,'uniformoutput',false)),results.fault.samples);
        [wRMS_min,i_max] = min(wRMSs);
        if results.settings.debug;fprintf('Mininmal rmsw: %f\n',wRMS_min);end
        fprintf('Minimal wRMS error: %f\n',wRMS_min);
        LeastSquaresFit = states{i_max} %#ok<NOPRT,NASGU>
        next_figure('Slip intensity')
        try
            subplot(2,1,1)
            intensity_plot(states,repeats,settings.T_min,settings.T_max,resolution_factor,plot_time);
            subplot(2,1,2)
            events_scatter_plot(states,repeats,fault,settings,plot_time,resolution_factor);
            drawnow
        catch e
            err(e)
        end
        next_figure('Measurements and simulations') 
        try
            simulations_plot(states,fault,repeats,i_max,confidence);
            drawnow
        catch e
            err(e)
        end
        next_figure('Parameters (1)')
        cols=3;rows=3;
        try
            subplot(rows,cols,1)
            histogram_plot(states,resolution_factor,repeats,@(x) (plot_time(2)-plot_time(1))/sum((x.times>plot_time(1))&(x.times<plot_time(2))), '$T_{mean}$',confidence);
            subplot(rows,cols,2)
            histogram_plot(states,resolution_factor,repeats,@(x) CV(x,plot_time),'CV',confidence)           
            subplot(rows,cols,3)
            histogram_plot(states,resolution_factor,repeats,@(x) x.rho_coll,'$\rho_{coll}$',confidence)
            subplot(rows,cols,4)
            histogram_plot(states,resolution_factor,repeats,@(x) x.rho,'$\rho$',confidence)
            subplot(rows,cols,5)
            histogram_plot(states,resolution_factor,repeats,@(x) sum((x.freq_switches>plot_time(1))&(x.freq_switches<plot_time(2))),'Number of switch points, $n$',confidence);
            subplot(rows,cols,6)
            histogram_plot(states,resolution_factor,repeats,@(x) sum((x.times>plot_time(1))&(x.times<plot_time(2))),'Number of earthquakes, $N$',confidence);
            subplot(rows,cols,7)
            histogram_plot(states,resolution_factor,repeats,@(x) x.times(end),'$T_{recent}$',confidence)
            subplot(rows,cols,8)
            histogram_plot(states,resolution_factor,repeats,@(x) x.T_init,'$T_{init}$',confidence)
            subplot(rows,cols,9)
            histogram_plot(states,resolution_factor,repeats,@(x) x.tau,'$\tau$',confidence)
            drawnow
        catch e
            err(e)
        end
        next_figure('Parameters (2)')
        cols=3;rows=3;
        try
            subplot(rows,cols,1)
            histogram_plot(states,resolution_factor,repeats, @(x) x.lambda_switches,'$\Lambda_{switches}$',confidence)
            subplot(rows,cols,2)
            histogram_plot(states,resolution_factor,repeats, @(x) x.freq_mean,'$f_{mean}$',confidence)
            subplot(rows,cols,3)
            histogram_plot(states,resolution_factor,repeats, @(x) x.freq_alpha,'$\alpha_{f}$',confidence)
            subplot(rows,cols,4)
            histogram_plot(states,resolution_factor,repeats,@(x) x.p_zero_freq,'$p_{f=0}$',confidence)
            subplot(rows,cols,5)
            histogram_plot(states,resolution_factor,repeats,@(x) x.Psi_sp,'$\Psi_{sp}$',confidence)
            subplot(rows,cols,6)
            histogram_plot(states,resolution_factor,repeats,@(x) x.Psi_mu,'$\Psi_{\mu}$',confidence)
            subplot(rows,cols,7)
            histogram_plot(states,resolution_factor,repeats,@(x) x.Lambda_sp,'$\Lambda_{sp}$',confidence)
            subplot(rows,cols,8)
            histogram_plot(states,resolution_factor,repeats,@(x) x.Lambda_mu,'$\Lambda_{mu}$',confidence)
            drawnow
        catch e
            err(e)
        end
        next_figure('Hazard assessment');
        try 
            next_earthquake_plot(states,repeats,settings);
            %correlation_plot(states,repeats,@(x) [x.parameters.correlation_length;x.parameters.rho;x.jumps(1);max(x.times)/1000;x.parameters.rho_coll;slip_between(x);x.T_init/1000;(plot_time(2)-plot_time(1))/1e3/sum((x.times>plot_time(1))&(x.times<plot_time(2)))],{'colength','modelerro','d_{init}','T_{recent}','\rho_{coll}','\Sigma_{d}','T_{init}','T_{recurrence}'});
            drawnow
        catch e
            err(e)
        end
        next_figure('Displacement history')
        try
            slip_history_plot(states,repeats,fault,settings,plot_time,i_max,confidence);
            drawnow
        catch e
            err(e)
        end
        next_figure('Recent earthquakes')
        try
            recent_earthquakes_plot(states,repeats,plot_time,'Time (kyr)',resolution_factor);
            drawnow
        catch e
            err(e)
        end
        next_figure('Convergence diagnostic')     
        try
            gelman_rubin_plot(results,settings,p_burn_in,plot_time);
            drawnow
        catch e
            err(e)
        end
    else
        for j=2:1+N_figs
            if ishandle(j)
                set(0,'CurrentFigure',j)
                clf
                text(0.0,0.5,'Not enough data','fontsize',50,'color','r')
                axis off
            end
        end
    end
    function clear_figs()
        for j=2:1+N_figs
            if ishandle(j)
                set(0,'CurrentFigure',j)
                clf
                text(-0.15,0.5,'Plot pending...','fontsize',50,'color','k')
                axis off
            end
        end
    end
    function next_figure(title)
       n_figure = n_figure+1;
       if ~ishandle(n_figure)
           figure(n_figure)
           clf
       else
           set(0,'CurrentFigure',n_figure);
           axis on
       end
       clf
       set(gcf,'Name',title,'NumberTitle','off')
    end
    function err(e)
        clf
        fprintf(2,'%s', getReport( e, 'extended', 'hyperlinks', 'on' ) );
        text(0.0,0.5,'Plot failed','fontsize',50,'color','r')
        axis off
    end
end
