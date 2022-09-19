clear all

Tvariable u [];
Tvariable d [];
Tvariable a [];

Tvariable Hess_ [6,6];

classname=class2minmaxCS('classname','tmp_minmax',...
                         'objective',(u+d)^2-a*d^2,...
                         'minOptimizationVariables',{u},...
                         'minConstraints',{u>-1,u<1},...
                         'maxOptimizationVariables',{d},...
                         'maxConstraints',{d>-1,d<1},...
                         'outputExpressions',{u,d,Hess_},...
                         'parameters',a,...
                         'scaleCost',0,...
                         'scaleInequalities',false,...
                         'solverVerboseLevel',4);

obj=feval(classname);

help obj

setP_a(obj,.5);

setV_u(obj,.1);
setV_d(obj,.1);

mu0=1;
maxIter=30;
saveIter=-1;
addEye2Hessian=[1e-9,1e-9,1e-9];
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);

out=getOutputs(obj)



