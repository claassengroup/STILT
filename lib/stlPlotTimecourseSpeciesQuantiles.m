function [ output_args ] = stlPlotTimecourseSpeciesQuantiles(sbmlModel, obs, opts, Q, iteration )
%STLPLOTTIMECOURSESPECIESQUANTILES Summary of this function goes here
%   Detailed explanation goes here

    hSQuantiles = 2;
    if (opts.ShowPlots)
        hFig=figure(hSQuantiles);
    else
        hFig=figure('Visible', 'on');
    end
    
    allTimes=unique(vertcat(obs.time{:}));
    N_cells = length(obs.cellNr);
    cols = lines(N_cells);
    qIdcs = [1; round(length(opts.LatentQuantiles)/2); length(opts.LatentQuantiles)];
    for sIdx = 1:length(opts.Species)
        subplot(length(opts.Species),1,sIdx);
        speciesName = opts.Species(sIdx).Name;
        for i=1:N_cells
            this_time = Q{i,2}; 
            plotTimeIdcs = 1:find(this_time,1,'last');
            if (~isempty(plotTimeIdcs))
                q1 = squeeze(Q{i,1}(plotTimeIdcs,sIdx,qIdcs(1)));
                q2 = squeeze(Q{i,1}(plotTimeIdcs,sIdx,qIdcs(2)));
                q3 = squeeze(Q{i,1}(plotTimeIdcs,sIdx,qIdcs(3)));
                shadedErrorBar(this_time(plotTimeIdcs), q2, [abs(q3 - q2)'; abs(q1-q2)'], {'Color', cols(i,:), 'markerfacecolor', cols(i,:), 'LineWidth', 1, 'LineStyle', '--'}, 1); % 'MarkerSize', 10, 'LineWidth', 2)
                xlim([allTimes(1) allTimes(end)]);
                hold on;
            end
            
            if (~isempty(opts.Species(sIdx).DataIdx))
                plot(obs.time{i}, obs.(speciesName){i}, '-s', 'markerfacecolor', cols(i,:), 'color', cols(i,:), 'markersize', 2)
            end
        end
        title(speciesName);
        xlabel('time')
        ylabel('abundance');
        hold off;
    end

    drawnow;
    if (opts.SavePlots)
        
        saveas(hFig, fullfile(opts.OutDir,'speciesquantiles', [int2str(iteration) '.png']));
    end
    if (~opts.ShowPlots)
        close(hFig);
    end

end