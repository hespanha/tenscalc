clear all

Tvariable u [];
Tvariable d [];
Tvariable a [];

Tvariable addEye2HessianU_ [];
Tvariable addEye2HessianD_ [];
Tvariable addEye2HessianEq_ [];

oExpr=struct('u',u,...
             'd',d,...
             'addEye2HessianU',addEye2HessianU_,...
             'addEye2HessianD',addEye2HessianD_,...
             'addEye2HessianEq',addEye2HessianEq_);

switch 4
  case 1
    % no constraints, diagonal Hessian
    oExpr.Lf_z=Tvariable('Lf_z_',[2]);
    oExpr.Lf_zz=Tvariable('Lf_zz_',[2,2]);
    oExpr.Hess=Tvariable('Hess_',[2,2]);

    objective=u^2-a*d^2;
    minConstraints={};
    maxConstraints={};

  case 2
    % no constraints, non-diagonal Hessian
    oExpr.Lf_z=Tvariable('Lf_z_',[2]);
    oExpr.Lf_zz=Tvariable('Lf_zz_',[2,2]);
    oExpr.Hess=Tvariable('Hess_',[2,2]);

    objective=(u+d+1)^2-a*(d-1)^2;
    minConstraints={};
    maxConstraints={};

  case 3
    oExpr.Hess=Tvariable('Hess_',[4,4]);

    objective=(u+d)^2-a*(d+2)^2;
    minConstraints={};
    maxConstraints={d>-1,d<1};

  case 4
    oExpr.Hess=Tvariable('Hess_',[4,4]);

    objective=(u+d)^2-a*d^2;
    minConstraints={u>-.5,u<.5};
    maxConstraints={};

  case 6
    oExpr.Hess=Tvariable('Hess_',[6,6]);

    objective=(u+d)^2-a*d^2;
    minConstraints={u>-1,u<1};
    maxConstraints={d>-1,d<1};
end

classname=class2minmaxCS('classname','tmp_minmax',...
                         'objective',objective,...
                         'minOptimizationVariables',{u},...
                         'minConstraints',minConstraints,...
                         'maxOptimizationVariables',{d},...
                         'maxConstraints',maxConstraints,...
                         'outputExpressions',oExpr,...
                         'parameters',a,...
                         'scaleCost',0,...
                         'addEye2Hessian',true,...
                         'adjustAddEye2Hessian',true,...
                         'scaleInequalities',false,...
                         'solverVerboseLevel',4);



obj=feval(classname);

%help obj

setP_a(obj,10);

setV_u(obj,.1);
setV_d(obj,.1);

mu0=1;
maxIter=30;
saveIter=-1;
addEye2Hessian=[1e-9,1e-9,1e-9];
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);

out=getOutputs(obj)
full(out.Hess)
eig(out.Hess)
full(out.Lf_zz)
