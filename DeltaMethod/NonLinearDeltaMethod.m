function results=NonLinearDeltaMethod(fun,results,varargin)
%%  Delta method for nonlinear relationship between estimated coefficients and parameters of interest
%   theta: estimated coefficients
%   phi: parameter of interest
%   relationship phi(theta)=phi

%%  DEFAULT PARAMETERS
pValues=[.9 .95 .99];
cellfrmt='%- 12.4f';
lbl=[];
lblstem=[];

%%  DYNAMIC READ IN
%   Parse input arguments
[opt, argList] = parseInputs(varargin, who);
%   Set options from the structure
setOptionsFromStruct(opt, who);
%   Set options from name-value pairs
setOptionsFromNameValuePairs(argList, who);

%%  PREP
%   read out estimates
if isfield(results,'coef')
    coef=results.coef;
    varcoef=results.varcoef;
elseif isfield(results,'betaHat')
    coef=results.betaHat;
    varcoef=results.betaVHat;
else
    error("structure does not have any known coefficients and covariance matrix.")
end
%   convert point estimates
phi=fun(coef);
%   derivative
Dfun=JacobianEst(@(p)fun(p),coef,'f0',phi);
varphi=Dfun*varcoef*Dfun';
stderrphi=diag(sqrt(varphi));

%%  TABLE 
%   check for variable names
if isempty(lblstem)
    lblstem="NLphi";
end
if isempty(lbl)
    lbl=lblstem+num2str((1:numel(phi))');
end
%   create table
SlimTable=CreateCoefficientTable(phi,stderrphi,'VarNames',lbl,'cellfrmt',cellfrmt,'CriticalValues',pValues);

%%  COLLECT OUTPUT
results.phi=phi;
results.varphi=varphi;
results.stderrphi=stderrphi;
if isfield(results,'tableStarsphi')
    results.tableStarsphi=[results.tableStarsphi SlimTable];
else
    results.tableStarsphi=SlimTable;
end
end