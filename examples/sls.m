% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.

% Many ways to skin a cat...
%
% This example is about finding the vector x that minimizes
%       J(x)=\|Ax-b\|^2
%

!rm -rf tmp* @tmp*
clear all

%generateSolver=@class2optimizeCS;
generateSolver=@cmex2optimizeCS;
compilerOptimization='-O0';

N=10000/25
n=800/25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construction of objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tvariable A [N,n];
Tvariable b N;
Tvariable x n;

y=A*x-b;

J=norm2(y)/N;

thisA=rand(N,n);
thisb=rand(N,1);
x0=.02*rand(n,1);
J0=norm(thisA*x0-thisb,2)^2/N;

fprintf('Cost for initial value = %g\n\n',J0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unconstrained minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A x-b\|^2
% w.r.t.   x

%% Using mldivide
fprintf('Unconstrained minimization using mldivide (\\)\n');
t00=clock;
xm=thisA\thisb;
fprintf('Elapsed time = %.3fms\n',1000*etime(clock,t00));

fprintf('\tcost = %g\n',norm(thisA*xm-thisb,2)^2/N);

fprintf('\nUnconstrained minimization\n');
classname=generateSolver('classname','tmp_minslsu',...
                         'objective',J,...
                         'optimizationVariables',{x},...
                         'outputExpressions',{J,x},...
                         'parameters',{A,b},...
                         'compilerOptimization',compilerOptimization,...
                         'solverVerboseLevel',3);
% Create object
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_b(obj,thisb);
% Initialize primal variables
setV_x(obj,x0);
% Solve optimization
mu0=1;
maxIter=30;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Justar,xustar]=getOutputs(obj);

%% Unconstrained minimization using a slack variable
% minimize v
% w.r.t.   x
% subject to v>=\|A x-b\|^2

Tvariable v [];
v0=J0+1;

fprintf('\nUnconstrained minimization using a slack variable\n');
classname=generateSolver('classname','tmp_minslsus1',...
                         'objective',v,...
                         'optimizationVariables',{x,v},...
                         'constraints',{v>=J},...
                         'outputExpressions',{J,x},...
                         'parameters',{A,b},...
                         'compilerOptimization',compilerOptimization,...
                         'solverVerboseLevel',3);
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_b(obj,thisb);
% Initialize primal variables
setV_x(obj,x0);
setV_v(obj,v0);
% Solve optimization
mu0=1;
maxIter=30;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Justar,xustar]=getOutputs(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A x-b\|^2
% w.r.t.   x
% subject to 0<=x<=.05
fprintf('\nConstrained minimization using ''primalDual'' method\n');
classname=generateSolver('classname','tmp_minslsc2',...
                         'method','primalDual',...
                         'objective',J,...
                         'optimizationVariables',{x},...
                         'constraints',{x>=0,x<=.05},...
                         'outputExpressions',{J,x},...
                         'parameters',{A,b},...
                         'compilerOptimization',compilerOptimization,...
                         'solverVerboseLevel',3);
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_b(obj,thisb);
% Initialize primal variables
setV_x(obj,x0);
% Solve optimization
mu0=1;
maxIter=30;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Jcstar,xcstar]=getOutputs(obj);

clf
plot([xustar,xcstar],'.')
legend('unconstrained','constrained')

%% Constrained optimization using quadprog
fprintf('\nConstrained minimization using quadprog\n');
options=optimset('algorithm','interior-point-convex');
t0=clock;
[xqpstar,fval,flag]=quadprog(thisA'*thisA/N,-thisA'*thisb/N,[],[],[],[],zeros(n,1),.05*ones(n,1),[],options);
fprintf('Elapsed time = %5.2fms (flag=%d, f-val =%.6f)\n',1000*etime(clock,t0),flag,norm(thisA*xqpstar-thisb,2)^2/N);
t0=clock;
[xqpstar,fval,flag]=quadprog(thisA'*thisA/N,-thisA'*thisb/N,[],[],[],[],zeros(n,1),.05*ones(n,1),[],options);
fprintf('Elapsed time = %5.2fms (flag=%d, f-val =%.6f)\n',1000*etime(clock,t0),flag,norm(thisA*xqpstar-thisb,2)^2/N);

%% Constrained optimization using CVX
if 1
    clear x
    cvx_begin
       variables x(n);
    
       minimize sum_square(thisA*x-thisb)/N;
       subject to 
          0<= x;
          x<= .05;
       t0=clock;
    cvx_end
    fprintf('Elapsed time = %5.2fs\n',etime(clock,t0));
end