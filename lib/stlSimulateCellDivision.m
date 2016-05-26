function [ X1, X2 ] = stlSimulateCellDivision( Xmother, opts )
% STLSIMULATCELLDIVISION copies/splits the content of mother cells to
% two daughter cells each.
%
% Xmother ... state of mother cell(s), cells-by-states
% opts ... STILT options structure
% X1,X2 ... new states of daugther cells, cells-by-states
%
% See also stlSetOptParam, stlSetOptSpecies

    BinomialApproxN = 100;

    X1 = zeros(size(Xmother));
    X2 = zeros(size(Xmother));

    copySpecies = strcmp({opts.Species.Division},'copy');
    X1(:,copySpecies) = Xmother(:,copySpecies);
    X2(:,copySpecies) = Xmother(:,copySpecies);

    binomialSpecies = strcmp({opts.Species.Division},'binomial');
    binomialPs = 0.5*binomialSpecies;
    for sIdx = find(binomialSpecies)
        if (~isempty(opts.Species(sIdx).DivisionParameters))
            strP = regexp(opts.Species(sIdx).DivisionParameters,'p=(\d+\.?\d*)', 'tokens');
            binomialPs(sIdx) = str2double(strP{1});
        end
    end
    binomialPs = repmat(binomialPs, [size(Xmother,1) 1]);
    binomialExact = repmat(binomialSpecies,[size(Xmother,1) 1]) & Xmother<BinomialApproxN;
    if (any(binomialExact(:)))
        sister1Values = binornd(Xmother(binomialExact),binomialPs(binomialExact));
        X1(binomialExact) = sister1Values;
        X2(binomialExact) = Xmother(binomialExact) - sister1Values;
    end
    binomialApprox = repmat(binomialSpecies,[size(Xmother,1) 1]) & Xmother>=BinomialApproxN;
    if (any(binomialApprox(:)))
        sister1Values = max(0, round(normrnd(Xmother(binomialApprox).*binomialPs(binomialApprox), sqrt(Xmother(binomialApprox)).*binomialPs(binomialApprox))));
        X1(binomialApprox) = sister1Values;
        X2(binomialApprox) = Xmother(binomialApprox) - sister1Values;
    end      
end

