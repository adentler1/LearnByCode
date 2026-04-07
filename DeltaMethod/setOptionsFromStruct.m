function setOptionsFromStruct(opt, varList)
%SETOPTIONSFROMSTRUCT Sets options from an options structure.
%
% Inputs:
%   opt     - Options structure from the input arguments.
%   varList - List of variable names in the main function's workspace.

if isempty(opt)
    return; % No options structure to process
end

for i = 1:length(varList)
    varName = varList{i};

    if isfield(opt, varName)
        newValue = opt.(varName);
        assignin('caller', varName, newValue);
    end
end
end
