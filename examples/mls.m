% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

% Many ways to skin a cat...
%
% This example is about finding the matrix X that minimizes
%       J(x)=\|A X-B\|^2
%

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

generateSolver=@class2optimizeCS;
%generateSolver=@cmex2optimizeCS;  % only for small problems (size/5)
compilerOptimization='-O0';

N=100;
n=30;
k=20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construction of objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tvariable A [N,n];
Tvariable B [N,k];

if 1
    Tvariable  X [n,k];
    Y=A*X-B;
else
    Tvariable  X [n*k];
    XX=reshape(X,[n,k]);
    Y=A*XX-B;
end

J=norm2(Y)/N;
%J=tprod(Y,[-1,-2],Y,[-1,-2]);
%J=sum((A*X-B).*(A*X-B),[1,2]);
%J=sum((A*X).*(A*X)+B.*B-2*(A*X).*B,[1,2]);

% initialize seed for repeatability
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

thisA=rand(N,n);
thisB=rand(N,k);
X0=.025+.02*rand(n,k);
J0=norm(thisA*X0-thisB,'fro')^2/N;

fprintf('Cost for initial value = %g\n\n',J0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unconstrained minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A X-B\|_F^2
% w.r.t.   X

fprintf('\nUnconstrained minimization using class2optimize\n');
classname=generateSolver('classname','tmp_minmlsu',...
                         'objective',J,...
                         'optimizationVariables',{X},...
                         'outputExpressions',{J,X},...
                         'parameters',{A,B},...
                         'compilerOptimization',compilerOptimization,...
                         'solverVerboseLevel',3);
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_B(obj,thisB);
% Initialize primal variables
setV_X(obj,X0);
% Solve optimization
mu0=1;
maxIter=100;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Justar,Xustar]=getOutputs(obj);

Tcalculus.clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A X-B\|_F^2
% w.r.t.   X
% subject to 0<=X<=.05

Tvariable A [N,n];
Tvariable B [N,k];

Tvariable  X [n,k];
Y=A*X-B;

J=norm2(Y)/N;

fprintf('\nConstrained minimization using class2optimize\n');
classname=generateSolver('classname','tmp_minmlsc',...
                         'objective',J,...
                         'optimizationVariables',{X},...
                         'constraints',{X>=0,X<=.05},...
                         'outputExpressions',{J,X},...
                         'parameters',{A,B},...
                         'compilerOptimization',compilerOptimization,...
                         'solverVerboseLevel',3);
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_B(obj,thisB);
% Initialize primal variables
setV_X(obj,X0);
% Solve optimization
mu0=1;
maxIter=100;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Jcstar,Xcstar]=getOutputs(obj);

clf
plot([Xustar(:),Xcstar(:)],'.')
legend('unconstrained','constrained')
