function [ s ] = cell2sym( c )

    allVars = symvar(strjoin(c,'+'));
    eval(['syms ' strjoin(allVars,' ') ';']);
    s = sym('S', size(c));
    for i = 1:numel(c)
        s(i) = sym(eval(c{i}));
    end
    
end

