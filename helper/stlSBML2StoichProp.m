function [ S, Sed, P ] = stlSBML2StoichProp( sbmlModel )
% STLSBML2STOICHPROP generates a stoichiometric matrix S and a cell array
% of propensitiy functions P from SBML model defintion

    SpeciesIDs = {sbmlModel.species.id};
    Ss = length(sbmlModel.species);
    Rs = length(sbmlModel.reaction);
    
    subsMap = containers.Map();
    for sIdx = 1:Ss
        subsMap(SpeciesIDs{sIdx}) = sprintf('X(:,%d)',sIdx);
    end
    for pIdx = 1:length(sbmlModel.parameter)
        subsMap(sbmlModel.parameter(pIdx).id) = sprintf('k(:,%d)',pIdx);
    end
    
    S = zeros(Ss, Rs);
    Sed = zeros(Ss, Rs);
    P = cell(Rs,1);
    for rIdx = 1:Rs
        reaction = sbmlModel.reaction(rIdx);
        [found,eductLoc]=ismember({reaction.reactant.species}, SpeciesIDs);
        if (any(~found))
            error('Unknown reactant in reaction %d', rIdx);
        end
        if (~isempty(eductLoc))
            S(eductLoc,rIdx) = S(eductLoc,rIdx)-[reaction.reactant.stoichiometry]';
            Sed(eductLoc,rIdx) = Sed(eductLoc,rIdx)-[reaction.reactant.stoichiometry]';
        end
        
        [found,productLoc]=ismember({reaction.product.species}, SpeciesIDs);
        if (any(~found))
            error('Unknown product in reaction %d', rIdx);
        end
        if (~isempty(productLoc))
            S(productLoc,rIdx) = S(productLoc,rIdx) + [reaction.product.stoichiometry]';
        end
        
        P{rIdx} = regexprep(vectorize(reaction.kineticLaw.formula), '([\w_]+)', '${stlSubsToken(subsMap, $1)}');
    end
end


