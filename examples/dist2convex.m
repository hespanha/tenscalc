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
% This example is about find the "convex combination" x that minimizes
%       J(x)= \|A x - b\|^2
% where sum(x)=1, x>=0

!rm -fr tmp* @tmp*
clear all

N=30;
d=9;

if 0
    % using equality contraint
    Tvariable lambda [N];

    optimizationVariables={lambda};
    constraints={sum(lambda,1)==1,lambda>=0};
else
    % using one less variable
    Tvariable lambda1 N-1;
    
    lambda=[lambda1;1-sum(lambda1,1)];
    optimizationVariables={lambda1};
    constraints={lambda>=0};
end

if 0
    % using original cost function
    Tvariable A [d,N];
    Tvariable b d;

    J=norm2(A*lambda-b);

    parameters={A,b};
else
    % expanding norm
    %   J=(lambda'*A'-b')*(A*lambda-b)
    %    =lambda'*(A'*A)*lambda-2*b'*A*lambda
    Tvariable AA  [N,N];
    Tvariable bA2 [N];

    J=lambda*AA*lambda-bA2*lambda;

    parameters={AA,bA2};
end


classname=cmex2optimizeCS('classname','tmp_mind2cvx',...
                          'objective',J,...
                          'optimizationVariables',optimizationVariables,...
                          'constraints',{sum(lambda,1)==1,lambda>=0},...
                          'outputExpressions',{J,lambda},...
                          'parameters',parameters,...
                          'solverVerboseLevel',3);

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
if exist('lambda1','var')
    lambda0=1/2/N*ones(N-1,1);
    setV_lambda1(obj,lambda0);
    lambda0=[lambda0;1-sum(lambda0)];
else
    lambda0=1/2/N*ones(N,1);
    setV_lambda(obj,lambda0);
end


% Solve optimization
mu0=1;
maxIter=30;
saveIter=-1;
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
[status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% Get outputs
[J_star,lambda_star]=getOutputs(obj);
J0=norm(b-A*lambda0,2)^2;
fprintf('>>> cost for initial value = %g\n',J0);
J=norm(b-A*lambda_star,2)^2;
fprintf('>>> optimal cost for unconstrained minimization = %g\n',J);
