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
% This example is about find the vector x that minimizes
%       J(x)=\|Ax-b\|^2
%

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

N=10000
n=800
m=40

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construction of objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tvariable A [N,n];
Tvariable b N;
Tvariable C [m,n];
Tvariable d m;
Tvariable x n;

eq=C*x-d;
y=A*x-b;

J=norm2(y);

thisA=rand(N,n);
thisb=rand(N,1);
thisC=rand(m,n);
thisd=rand(m,1);
x0=.01*rand(n,1);

s=norm(thisb);
thisb=thisb/s;
thisA=thisA/s;

s=norm(thisd);
thisd=thisd/s;
thisC=thisC/s;

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

fprintf('\tcost = %g\n',norm(thisA*xm-thisb,2)^2);

fprintf('\nUnconstrained minimization\n');
classname=class2optimizeCS('classname','tmp_minslsu',...
                           'objective',J,...
                           'optimizationVariables',{x},...
                           'parameters',{A,b},...
                           'outputExpressions',{J,x},...
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
maxIter=100;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Justar,xustar]=getOutputs(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained minimization with linear equality constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A x-b\|^2
% w.r.t.   x
% subject to C x == d

%% Using class2optimize
fprintf('\nMinimization with linear equality constraints\n');
classname=class2optimizeCS('classname','tmp_minslsceq',...
                           'objective',J,...
                           'optimizationVariables',{x},...
                           'constraints',{eq==0},...
                           'parameters',{A,b,C,d},...
                           'outputExpressions',{J,x},...
                           'solverVerboseLevel',3);
% Create object
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_b(obj,thisb);
setP_C(obj,thisC);
setP_d(obj,thisd);
% Initialize primal variables
setV_x(obj,x0);
% Solve optimization
mu0=1;
maxIter=100;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Justar,xustar]=getOutputs(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained minimization with inequality contraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A x-b\|^2
% w.r.t.   x
% subject to -.02<=x<=.02

fprintf('\nConstrained minimization with inequality contraints\n');
classname=class2optimizeCS('classname','tmp_minslscin',...
                           'objective',J,...
                           'optimizationVariables',{x},...
                           'constraints',{x>=-.02,x<=.02},...
                           'parameters',{A,b,C,d},...
                           'outputExpressions',{J,x},...
                           'solverVerboseLevel',3);
% Create object
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_b(obj,thisb);
setP_C(obj,thisC);
setP_d(obj,thisd);
% Initialize primal variables
setV_x(obj,x0);
% Solve optimization
mu0=1;
maxIter=100;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Justar,xcstar]=getOutputs(obj);

min(xcstar)
max(xcstar)
norm(thisC*xcstar-thisd)


subplot(2,1,1)
plot([xustar,xcstar],'.')
legend('unconstrained','constrained inequality')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained minimization with equality & inequality contraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A x-b\|^2
% w.r.t.   x
% subject to -.02<=x<=.02
% subject to C x == d

fprintf('\nConstrained minimization with equality & inequality contraints\n');
classname=class2optimizeCS('classname','tmp_minslsceqin',...
                           'objective',J,...
                           'optimizationVariables',{x},...
                           'constraints',{eq==0,x>=-.02,x<=.02},...
                           'parameters',{A,b,C,d},...
                           'outputExpressions',{J,x},...
                           'solverVerboseLevel',3);
% Create object
obj=feval(classname);
% Set parameters
setP_A(obj,thisA);
setP_b(obj,thisb);
setP_C(obj,thisC);
setP_d(obj,thisd);
% Initialize primal variables
setV_x(obj,x0);
% Solve optimization
mu0=1;
maxIter=30;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[Justar,xcstar]=getOutputs(obj);

min(xcstar)
max(xcstar)
norm(thisC*xcstar-thisd)

subplot(2,1,2)
plot([xustar,xcstar],'.')
legend('unconstrained','constrained equality & inequality')


