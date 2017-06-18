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
% This example is about finding the matrix X that minimizes
%       J(x)=\|A X-B\|^2
%

!rm -rf tmp* @tmp*
clear all
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

generateSolver=@class2optimizeCS;
%generateSolver=@cmex2optimizeCS;  % only for small problems (/5)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize \|A X-B\|_F^2
% w.r.t.   X
% subject to 0<=X<=.05

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
