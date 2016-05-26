function [ hFig ] = stlPlotSampleFrequencies( samples, sbmlModel, opts, iteration )
% STLPLOTSAMPLEFREQUENCIES plots a histgram of the re-sampling frequencies
% of individual particles
% 
% See also stlPlotPosterior, stlPlotPosterior

    hFrequency = 3;
    if (opts.ShowPlots)
        hFig=figure(hFrequency);
        clf
    else
        hFig=figure('Visible', 'off');
    end

    freqs = tabulate(samples);
    freqs2 = tabulate(freqs(:,2));
    stem(freqs2(:,1),freqs2(:,2)/opts.ParticleNr);
    xlabel('Sampling frequency');
    ylabel('Fraction of particles');
    title('Particle sampling frequencies');
    
    drawnow;
    if (opts.SavePlots)
        saveas(hFig, fullfile(opts.OutDir, 'frequencies', [int2str(iteration) '.png']));
    end
    if (~opts.ShowPlots)
        close(hFig);
    end 
end

