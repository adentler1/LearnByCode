function [opt, argList] = parseInputs(inputArgs, varList)
%PARSEINPUTS Parses input arguments for structures and name-value pairs.
%
% Inputs:
%   inputArgs - Cell array of varargin from the main function.
%   varList   - List of variable names in the main function's workspace.
%
% Outputs:
%   opt     - Options structure if provided in the input arguments.
%   argList - Cell array containing name-value pairs.

% Initialize outputs
opt = [];
argList = [];

% Check if there are any additional arguments
if isempty(inputArgs)
    return; % No additional arguments to parse
end

% Check if the first argument is a structure (options)
if isstruct(inputArgs{1})
    opt = inputArgs{1};
    inputArgs(1) = []; % Remove the processed argument
end

% Remaining arguments should be name-value pairs
if mod(length(inputArgs), 2) ~= 0
    error('Expected name-value pairs after the options structure.');
end

% Validate and collect name-value pairs
for i = 1:2:length(inputArgs)
    name = inputArgs{i};
    value = inputArgs{i + 1};

    % Validate name
    if ~ischar(name) && ~isstring(name)
        error('Name must be a character vector or string scalar.');
    end

    % Check if name is a valid variable in the main function
    if ~ismember(name, varList)
        warning(['Ignoring unknown option: ', char(name)]);
        continue;
    end

    % Add to argList
    argList = [argList, {name, value}];
end
end
