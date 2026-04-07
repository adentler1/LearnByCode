function setOptionsFromNameValuePairs(argList, varList)
%SETOPTIONSFROMNAMEVALUEPAIRS Updates variables with name-value pairs.
%
% Inputs:
%   argList - Cell array of name-value pairs.
%   varList - List of variable names in the main function's workspace.

for i = 1:2:length(argList)
    name = argList{i};
    value = argList{i + 1};

    if ismember(name, varList)
        assignin('caller', name, value);
    else
        warning(['Ignoring unknown option: ', char(name)]);
    end
end
end
