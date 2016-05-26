function [ hFig ] = stlPlotPosterior( c, sbmlModel, opts, iteration )
% STLPLOTPOSTERIOR plots histograms of posterior P(Parameter|Data) per parameter 
% defined by particles c (particles-by-parameters) 
% 
% See also stlPlotTracking, stlPlotSampleFrequencies

    hPosterior = 1;
    if (opts.ShowPlots)
        hFig=figure(hPosterior);
        clf
    else
        hFig=figure('Visible', 'off');
    end
    
    % posterior for parameters
    for pIdx=1:length(sbmlModel.parameter)
        subplot(opts.SubplotRows,opts.SubplotColumns,pIdx); 
        hold on;
        [y1,x]=hist(c(:,pIdx),200);
        y1=y1/sum(y1);
        bar(x,y1); 
        prior = opts.Parameters(pIdx).GammaPrior;
        mu=gamstat(prior(1), 1/prior(2));
        h=line([mu,mu], get(gca,'YLim'));      
        set(h,'Color','r','LineStyle','--')
        if (opts.TrueParamValuesInModel)
            h=line(sbmlModel.parameter(pIdx).value * [1 1], get(gca,'YLim'));
            set(h,'Color','b','LineWidth',2)
        end
        [x,y2]=fplot( @(x) gampdf(x, prior(1), 1/prior(2)), get(gca,'XLim'));
        plot(x,y2/sum(y2), 'r--', 'LineWidth', 2)
        title(sbmlModel.parameter(pIdx).name, 'FontSize', 24)
    end
    drawnow;
    if (opts.SavePlots)
        if (~isempty(iteration))
            sIteration = int2str(iteration);
        else
            sIteration = 'posterior';
        end
        saveas(hFig, fullfile(opts.OutDir,'posterior', [sIteration '.png']));
    end
    if (~opts.ShowPlots)
        close(hFig);
    end
end

