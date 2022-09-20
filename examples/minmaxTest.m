clear all

Tvariable u [];
Tvariable d [];
Tvariable x [];
Tvariable a [];

Tvariable addEye2HessianU_ [];
Tvariable addEye2HessianD_ [];
Tvariable addEye2HessianEq_ [];

minOptimizationVariables={u};
maxOptimizationVariables={d};

oExpr=struct('u',u,...
             'd',d,...
             'addEye2HessianU',addEye2HessianU_,...
             'addEye2HessianD',addEye2HessianD_,...
             'addEye2HessianEq',addEye2HessianEq_);

switch 5.5
  case 1
    % no constraints, diagonal Hessian
    oExpr.Hess=Tvariable('Hess_',[2,2]);
    oExpr.HessD=Tvariable('HessD_',[1,1]);
    oExpr.HessU=Tvariable('HessU_',[2,2]);

    objective=u^2-2*d^2;
    minConstraints={};
    maxConstraints={};

  case 2
    % no constraints, non-diagonal Hessian
    oExpr.Hess=Tvariable('Hess_',[2,2]);
    oExpr.HessD=Tvariable('HessD_',[1,1]);
    oExpr.HessU=Tvariable('HessU_',[2,2]);

    objective=(u+d+1)^2-2*(d-1)^2;
    minConstraints={};
    maxConstraints={};

  case 2.5
    % equality constraint, non-diagonal Hessian
    oExpr.Hess=Tvariable('Hess_',[4,4]);
    oExpr.HessD=Tvariable('HessD_',[3,3]);
    oExpr.HessU=Tvariable('HessU_',[4,4]);

    maxOptimizationVariables{end+1}=x;

    objective=(x+1)^2-2*(d-1)^2;
    minConstraints={};
    maxConstraints={x==u+d};

  case 3
    oExpr.Hess=Tvariable('Hess_',[4,4]);
    oExpr.HessD=Tvariable('HessD_',[3,3]);
    oExpr.HessU=Tvariable('HessU_',[4,4]);
    
    objective=(u+d)^2-2*(d+2)^2;
    minConstraints={};
    maxConstraints={d>-1,d<1};

  case 3.5
    oExpr.Hess=Tvariable('Hess_',[6,6]);
    oExpr.HessD=Tvariable('HessD_',[5,5]);
    oExpr.HessU=Tvariable('HessU_',[6,6]);
    
    maxOptimizationVariables{end+1}=x;

    objective=(x)^2-2*(d+2)^2;
    minConstraints={};
    maxConstraints={d>-1,d<1,x==u+d};

  case 4
    oExpr.Hess=Tvariable('Hess_',[4,4]);
    oExpr.HessD=Tvariable('HessD_',[1,1]);
    oExpr.HessU=Tvariable('HessU_',[4,4]);

    objective=(u+d+1)^2-2*d^2;
    minConstraints={u>-.25,u<.25};
    maxConstraints={};

  case 5
    oExpr.Hess=Tvariable('Hess_',[6,6]);
    oExpr.HessD=Tvariable('HessD_',[3,3]);
    oExpr.HessU=Tvariable('HessU_',[6,6]);

    objective=(u+d)^2-2*(d+2)^2;
    minConstraints={u>-2,u<2};
    maxConstraints={d>-1,d<1};
  
  case 5.5
    oExpr.Hess=Tvariable('Hess_',[8,8]);
    oExpr.HessD=Tvariable('HessD_',[5,5]);
    oExpr.HessU=Tvariable('HessU_',[8,8]);

    maxOptimizationVariables{end+1}=x;

    objective=(x)^2-2*(d+2)^2;
    minConstraints={u>-2,u<2};
    maxConstraints={d>-1,d<1,x==u+d};
end

oExpr.J=objective;
classname=class2minmaxCS('classname','tmp_minmax',...
                         'objective',objective,...
                         'minOptimizationVariables',minOptimizationVariables,...
                         'minConstraints',minConstraints,...
                         'maxOptimizationVariables',maxOptimizationVariables,...
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

setP_a(obj,3);

setV_u(obj,.1);
setV_d(obj,.1);
if numel(maxOptimizationVariables)>1
    setV_x(obj,.1);
end

mu0=1;
maxIter=30;
saveIter=-1;
addEye2Hessian=[1e-9,1e-9,1e-9];
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter),addEye2Hessian);

out=getOutputs(obj)
full(out.Hess)
full(out.HessU)
full(out.HessD)
eig(out.Hess)
