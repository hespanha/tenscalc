% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
!rm -fr @tmp* tmp*;

%% Unconstrained

Tvariable P [3,3];
Tvariable x [3];

J=x*(P*x);

out=Tvars2optimizeCS('objective',J,...
                     'optimizationVariables',{x},...
                     'sensitivityVariables',{x},...
                     'constraints',{})
oExpr={J,x,full(out.Hess),full(out.Du1),full(out.DfDu1),full(out.D2fDu1)};

classname=cmex2optimizeCS('classname','tmp_test',...
                          'objective',J,...
                          'optimizationVariables',{x},...
                          'sensitivityVariables',{x},...
                          'constraints',{},...
                          'parameters',{P},...
                          'outputExpressions',oExpr,...
                          'addEye2Hessian',1e-10,...
                          'umfpack',false,...
                          'solverVerboseLevel',3);
obj=feval(classname);

P=rand(3);
P=P'*P;
setP_P(obj,P);
setV_x(obj,rand(3,1))

mu0=1;
maxIter=200;
saveIter=0;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));

[oExpr{:}]=getOutputs(obj);

oExpr{:},2*P,


clear all;

%% Constrained

Tvariable P [3,3];
Tvariable x [3];
Tvariable z [3];

J=x*(P*x);

out=Tvars2optimizeCS('objective',J,...
                     'optimizationVariables',{x,z},...
                     'sensitivityVariables',{z},...
                     'constraints',{x==z,x>=1})
oExpr={J,x,z,full(out.Hess),full(out.Du1),full(out.DfDu1),full(out.D2fDu1)};

classname=class2optimizeCS('classname','tmp_test',...
                          'objective',J,...
                          'optimizationVariables',{x,z},...
                          'sensitivityVariables',{z},...
                          'constraints',{x==z,z<=10},...
                          'parameters',{P},...
                          'outputExpressions',oExpr,...
                          'addEye2Hessian',1e-10,...
                          'umfpack',false,...
                          'solverVerboseLevel',3);
obj=feval(classname);

P=rand(3);
P=P'*P;
setP_P(obj,P);
setV_x(obj,5+rand(3,1))
setV_z(obj,5+rand(3,1))

mu0=1;
maxIter=200;
saveIter=0;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));

[oExpr{:}]=getOutputs(obj);

oExpr{:},2*P,

