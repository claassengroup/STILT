function [ replacement ] = stlSubsToken(subsMap, token)
    if (isKey(subsMap, token))
        replacement = subsMap(token);
    else
        replacement = token;
    end
end

