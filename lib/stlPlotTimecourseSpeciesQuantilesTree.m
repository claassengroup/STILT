function stlPlotTimecourseSpeciesQuantilesTree(sbmlModel, obs, opts, Q, iteration )
%STLPLOTTIMECOURSESPECIESQUANTILES Summary of this function goes here
%   Detailed explanation goes here

    N_cells = length(obs.cellNr);
    cellHeight = 300;
    genWidth = 300;
    generations = ceil(log2(N_cells+1));
    maxVertCells = 2^(generations-1);

    timepoints = sort(unique(vertcat(obs.time{:})));
    maxTimepoint = max(timepoints);
    cols = lines(N_cells);
    allQs = vertcat(Q{:,1});
    axisHeight = 1/maxVertCells*0.9;
    
    
    for sIdx = 1:length(opts.Species)

        speciesName = opts.Species(sIdx).Name;
        
        if (~isempty(opts.Species(sIdx).PlotQuantiles))
            [~,qIdcs] = ismember(opts.Species(sIdx).PlotQuantiles, opts.LatentQuantiles);
        else
            qIdcs = unique([1; round(length(opts.LatentQuantiles)/2); length(opts.LatentQuantiles)]);
        end
        maxTopQs = squeeze(max(allQs(:,:,qIdcs(end)),[],1));

        hFig=figure('position', [1 1 generations*genWidth maxVertCells*cellHeight], 'Color', 'white','visible','on');
    
        genT = repmat([Inf, -Inf],[generations 1]);
        for i=1:N_cells
            generation = ceil(log2(i+1));
            curMin = min(Q{i,2});
            curMax = max(Q{i,2});
            if (curMin < genT(generation,1))
                genT(generation,1) = curMin;
            end
            if (curMax > genT(generation,2))
                genT(generation,2) = curMax;
            end            
        end
        
        for i=1:N_cells
            generation = ceil(log2(i+1));
            curGenerationSize = 2^(generation-1);
            withinGenerationIdx = i-2^(generation-1)+1;
            genStartIdx = 2^generations/(2*curGenerationSize)-1;
            genStepIdx = 2^generations/curGenerationSize;
            
            this_time = Q{i,2}; 
            axes('outerposition', [this_time(1)/maxTimepoint (genStartIdx+(withinGenerationIdx-1)*genStepIdx)/2^generations diff(genT(generation,:))/maxTimepoint axisHeight]);

            if (length(qIdcs) == 3)
                q1 = squeeze(Q{i,1}(:,sIdx,qIdcs(1)));
                q2 = squeeze(Q{i,1}(:,sIdx,qIdcs(2)));
                q3 = squeeze(Q{i,1}(:,sIdx,qIdcs(3)));

                shadedErrorBar(this_time, q2, [abs(q3 - q2)'; abs(q1-q2)'], {'Color', cols(i,:), 'markerfacecolor', cols(i,:), 'LineWidth', 1, 'LineStyle', '--'}, 1); % 'MarkerSize', 10, 'LineWidth', 2)
            else
                plot(this_time, squeeze(Q{i,1}(:,sIdx,qIdcs)),'Color', cols(i,:), 'markerfacecolor', cols(i,:), 'LineWidth', 1, 'LineStyle', '-');
            end
            hold on;
            if (~isempty(opts.Species(sIdx).DataIdx))
                plot(obs.time{i}, obs.(speciesName){i}, '-s', 'markerfacecolor', cols(i,:), 'color', cols(i,:), 'markersize', 2)
            end
            xlim(genT(generation,:));
            ylim([0 max(maxTopQs(sIdx),1)]);
            xlabel('time')
            ylabel('abundance');
            title(sprintf('Cell %d', i));
        end
        suptitle(speciesName);
        if (opts.SavePlots)
            if ~isempty(iteration)
                outFile = [speciesName '_' int2str(iteration) '.png'];
            else
                outFile = [speciesName '.png'];
            end
            filePath = fullfile(opts.OutDir,'speciesquantilestree', outFile);
            %saveas(hFig, filePath);
            hFig.PaperPositionMode = 'auto';
            print(filePath,'-dpng','-r150')            
        end
        if (opts.ShowPlots)
            set(hFig, 'visible', 'on');
        else
            close(hFig);
        end
    end
end