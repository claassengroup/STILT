function [ hFig ] = stlPlotTracking( X, obs, cells, samples, data_times, observedSpecies, sbmlModel, opts, iteration )
% STLPLOTTRACKING plots observed and predicted densities 
% 
% See also stlPlotPosterior, stlPlotSampleFrequencies

    hTracking = 4;
    if (opts.ShowPlots)
        hFig=figure(hTracking);
    else
        hFig=figure('Visible', 'off');
    end
    hold on;

    dts = diff(data_times);
    dt = dts(iteration);
    cols = jet(length(cells));
    for sIdx = 1:length(observedSpecies)
        subplot(ceil(sqrt(length(observedSpecies))), floor(sqrt(length(observedSpecies))), sIdx);
        speciesIdx=observedSpecies(sIdx);
        speciesName=opts.Species(observedSpecies(sIdx)).Name;
        for k=1:length(cells)
            cell_id = find(obs.cellNr==cells(k));
            ix = find(obs.time{cell_id} == data_times(iteration+1));
            plot(dt*iteration, obs.(speciesName){cell_id}(ix), '*', 'MarkerSize', 10, 'Color', cols(k,:));  hold on;
            [f,x] = ksdensity(X{k}(:,speciesIdx));
            plot(dt*iteration + f/sum(f), x, 'Color', cols(k,:))
            [f,x] = ksdensity(X{k}(samples,speciesIdx));
            plot(dt*iteration + f/sum(f), x, 'k')                      
        end
        xlabel('time');
        ylabel('abundance');
        title(opts.Species(speciesIdx).Name);
    end
    drawnow;
    if (opts.SavePlots)
        saveas(hFig, fullfile(opts.OutDir,'tracking', [int2str(iteration) '.png']));
    end
    if (~opts.ShowPlots)
        close(hFig);
    end    
end

