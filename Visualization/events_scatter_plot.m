function events_scatter_plot(states,repeats,fault,settings,plot_time,resolution_factor)
    Z = cellfun(@(state) [state.times(2:end)/1000.0,state.jumps(2:end)],states,'UniformOutput',false);
    Z = cell2mat(Z);
    X = Z(:,1);
    Y = Z(:,2);
    weightedhist(X,Y,repeats,fault,settings,resolution_factor)
    ylabel('Displacement size (cm)','interpreter','latex')
    xlabel('Time (kyr)','interpreter','latex')
    axis([plot_time(1)/1000.0,plot_time(2)/1000.0,0,1.2*max(Y)]);
    hold on;
    if isfield(fault.truth,'times') && isstruct(fault.truth)
        plot(fault.truth.times/1000.0,fault.truth.jumps,'b.','MarkerSize',15);
        legend({'Truth'})
    end
end

function weightedhist(X,Y,W,fault,settings,res_factor)%https://stackoverflow.com/questions/6777609/fast-2dimensional-histograming-in-matlab/6778256#6778256
    xbins = linspace(settings.T_min/1000.0,settings.T_max/1000.0,8*res_factor*ceil(sum(W)^(1/4)));
    ybins = linspace(0,fault.parameters.H_sc,2*res_factor*ceil(sum(W)^(1/4)));
    xNumBins = numel(xbins); 
    yNumBins = numel(ybins);
    Xi = round( interp1(xbins, 1:xNumBins, X, 'linear', 'extrap') );
    Yi = round( interp1(ybins, 1:yNumBins, Y, 'linear', 'extrap') );
    Xi = max( min(Xi,xNumBins), 1);
    Yi = max( min(Yi,yNumBins), 1);
    H = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
    imagesc(xbins, ybins, H)
    axis on
    set(gca,'YDir','normal')
    colormap(flipud(gray));
end
