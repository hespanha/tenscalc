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
% This example is about finding the "convex combination" x that minimizes
%       J(x)= \|A x - b\|^2
% where sum(x)=1, x>=0

clear all
!rm -rf toremove.m tmp* @tmp*;

solverGeneration=@cmex2optimizeCS;

N=100;
d=9;

if 1
    % using equality contraint
    Tvariable x [N];

    optimizationVariables={x};
    constraints={sum(x,1)==1,x>=0};
else
    % using one less variable
    Tvariable x1 N-1;
    
    x=[x1;1-sum(x1,1)];
    optimizationVariables={x1};
    constraints={x>=0};
end

if 0
    % using original cost function
    Tvariable A [d,N];
    Tvariable b d;

    J=norm2(A*x-b);

    parameters={A,b};
else
    % expanding norm
    %   J=(x'*A'-b')*(A*x-b)
    %    =x'*(A'*A)*x-2*b'*A*x
    Tvariable AA  [N,N];
    Tvariable bA2 [N];

    J=x*AA*x-bA2*x;

    parameters={AA,bA2};
end


[classname,code]=solverGeneration('classname','tmp_mind2cvx',...
                                  'objective',J,...
                                  'optimizationVariables',optimizationVariables,...
                                  'constraints',constraints,...
                                  'outputExpressions',{J,x},...
                                  'parameters',parameters,...
                                  'compilerOptimization','-O1',...
                                  'profiling',true,...
                                  'solverVerboseLevel',2);

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

A=randn(d,N);
b=randn(d,1);

obj=feval(classname);

% Set parameters
if exist('AA','var')
    setP_AA(obj,A'*A);
    setP_bA2(obj,2*A'*b);
else
    setP_A(obj,A);
    setP_b(obj,b);
end
% Initialize primal variables
if exist('x1','var')
    x0=1/2/N*ones(N-1,1);
    setV_x1(obj,x0);
    x0=[x0;1-sum(x0)];
else
    x0=1/2/N*ones(N,1);
    setV_x(obj,x0);
end


% Solve optimization
mu0=1;
maxIter=30;
saveIter=-1;
clear iter time;
nRepeats=1000;
for i=1:nRepeats
    [status,iter(i),time(i)]=solve(obj,mu0,int32(maxIter),int32(saveIter));
end
fprintf('iter:   mean=%.2f, min=%.2f, max=%.2f\n',...
        mean(iter),min(iter),max(iter));
fprintf('time [us]:   median=%.2f, min=%.2f, max=%.2f\n',...
        1e6*median(time),1e6*min(time),1e6*max(time));
