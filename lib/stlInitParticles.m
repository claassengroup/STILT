function [ X ] = stlInitParticles( P, X0, S0, opts )
% STLINITPARTICLES initialized the particles' latent states from the
% Matlab expression defined in options.Species.Init. 
%
% Variables available for the expression:
% <ObsSpeciesName>1 ... First measured value for observed species
% <ObsSpeciesName>Sigma1 ... First measurement error of observed species
% E.g. Observed 'Protein' initialized 'normrnd(Protein1, ProteinSigma1, [P, 1])'
%
% Initialization expressions are evaluated in the species' order of the options
% structure and added to the environment as
% <SpeciesName> ...  P x 1 vector of values of beforehand assigned species
% E.g. if species 'DNA_on' was defined before, 'DNA_off' can be set by '1-DNA_on'
%
% See also stlOptions

    observedSpecies = find(cellfun(@(c) c==true, {opts.Species.Observed}));
    
    assert(length(observedSpecies)==length(X0), 'Number of observed species ~= length of datapoints');
    for sIIdx = 1:length(observedSpecies)
        sIdx = observedSpecies(sIIdx);
        speciesID = opts.Species(sIdx).ID;
        eval(sprintf('%s1=%d;',speciesID,X0(sIIdx)));
        eval(sprintf('%sSigma1=%d;',speciesID,S0(sIIdx)));
    end
    
    clear sIIdx observedSpecies;
    
    X = zeros(P, length(opts.Species));
    for sIdx = 1:length(opts.Species)
        speciesID = opts.Species(sIdx).ID;
        expression = opts.Species(sIdx).Init;
        eval(sprintf('%s=%s;',speciesID,expression));
        X(:,sIdx) = eval(sprintf('%s',speciesID));
    end
    
end

