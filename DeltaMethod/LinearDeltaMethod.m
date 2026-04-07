function results=LinearDeltaMethod(R,r,results)
%%  Delta method for linear relationship between estimated coefficients and parameters of interest
%   theta: estimated coefficients
%   phi: parameter of interest
%   relationship phi(theta)=phi
%   linear relationship: R*theta-r=phi

%%  identification
if size(R,2)~=results.k
    error('R bad!')
end
results.phi=R*results.coef-r;
results.varphi=R*results.varcoef*R';
results.stderrphi=diag(sqrt(results.varphi));

%%  TABLE 
%   table of variables with significance measure
lbl="Lphi"+num2str((1:size(R,1))');
tab=CreateCoefficientTable(results.phi,results.stderrphi, ...
    'VarNames',lbl,'cellfrmt','%- 12.4f','CriticalValues',results.pValues);
%%  COLLECT OUTPUT
if isfield(results,'tableStarsphi') 
results.tableStarsphi=[results.tableStarsphi tab];
else
results.tableStarsphi=tab;
end
end
