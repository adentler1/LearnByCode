function [SlimTable,FullTable]=CreateCoefficientTable(PE,SD,varargin)
%CREATECOEFFICIENTTABLE  Formats point estimates with significance stars.
%
%   [SlimTable, FullTable] = CreateCoefficientTable(PE, SD, ...)
%
%   Inputs:
%       PE  - Column vector of point estimates
%       SD  - Column vector of standard deviations (standard errors)
%
%   Optional name-value pairs:
%       'VarNames'       - String array of variable names (default: "beta1", ...)
%       'CriticalValues' - Significance levels (default: [.9 .95 .99])
%       'cellfrmt'       - Number format string (default: '%- 12.4f')
%       'Sidedness'      - "one-sided" (default) or "two-sided"
%       'TStats'         - Pre-computed t-statistics (default: PE./SD)
%       'PValues'        - Pre-computed p-values (default: from t-stats)
%
%   Outputs:
%       SlimTable - 1-row table with formatted estimate strings (stars + SE)
%       FullTable - k-row table with PE, SD, t-stats, p-values, formatted strings

PValues=[];
TStats=[];
CriticalValues=[.9 .95 .99];
cellfrmt='%- 12.4f';
VarNames=[];
ColNames=["PointEstimates","StandardDeviation","T-Statistic","PValues","TableValues"];
Sidedness="one-sided"; %"two-sided"

%%  DYNAMIC READ IN
%   Parse input arguments
[opt, argList] = parseInputs(varargin, who);
%   Set options from the structure
setOptionsFromStruct(opt, who);
%   Set options from name-value pairs
setOptionsFromNameValuePairs(argList, who);

%%  Prep
%   normalize the critical values
if any(CriticalValues>1)
    error("CriticalValues are between 0 and 1.")
end
CriticalValues(CriticalValues>.5)=1-CriticalValues(CriticalValues>.5);
CriticalValues=CriticalValues(:)';
%   number of coefficients
k=numel(PE);
%   calculate the implied p-value from a one-sided t statistic with an
%   infinite degree of freedom
if isempty(TStats)
    TStats=NaN(k,1);
    for i1=1:k
        TStats(i1)=PE(i1)/SD(i1);
    end
end
if isempty(PValues)
    PValues=NaN(k,1);
    for i1=1:k
        PValues(i1)=tcdf(abs(TStats(i1)),inf,'upper');
    end
    if Sidedness=="two-sided"
        PValues=PValues*2;
    end
end
%   number of stars implies
Stars=NaN(k,1);
for i1=1:k
    Stars(i1)=sum(PValues(i1)<=CriticalValues);
end
%   commented lines
C=strings(k,1);
for i1=1:k
    if Stars(i1)>0
        tmp=string(repmat('*',1,Stars(i1)));
    else
        tmp="";
    end
    if isfinite(SD(i1))
    C(i1)=num2str(PE(i1),cellfrmt)+tmp+newline+"("+ num2str(SD(i1),cellfrmt)+")";
    else
    C(i1)=num2str(PE(i1),cellfrmt)+tmp+newline+"(p="+ num2str(PValues(i1),cellfrmt)+")";

    end
end
%   variable names
if isempty(VarNames)
    VarNames=strings(k,1);
    for i1=1:k
        VarNames(i1)="beta"+num2str(i1);
    end
end

%%  table of variables with significance measure
FullTable=[table(PE(:),SD(:),TStats(:),PValues(:),C,...
    'VariableNames',ColNames,...
    'RowNames',VarNames)....
     ];
SlimTable=array2table(C','VariableNames',VarNames);
end
