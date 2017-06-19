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

clear all
!rm -rf toremove.m tmp* @tmp*;

solverGeneration=@cmex2optimizeCS;

N=30;
d=9;

if 1
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


classname=solverGeneration('classname','tmp_mind2cvx',...
                           'objective',J,...
                           'optimizationVariables',optimizationVariables,...
                           'constraints',constraints,...
                           'outputExpressions',{J,lambda},...
                           'parameters',parameters,...
                           'compilerOptimization','-O1',...
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

% >> [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% ipmPD_CSsolver.c (skipAffine=0,delta=3,allowSave=0): 300 primal variables, 1 eq. constr., 300 ineq. constr.
% Iter   cost      |grad|      |eq|    inequal     dual      gap       mu      alphaA    sigma     alphaS   time [us]
%  30:<-mx tol->   1.00e-04  1.00e-04                      1.00e-05  1.67e-08
%   1:  1.820e-01  6.13e+02  5.00e-01  1.67e-03  6.00e+02  3.00e+02  1.00e+00  4.94e-01  -sigma-   1.00e+00 83282.2us
%   2:  2.753e-01  4.53e-12  1.20e-06  3.14e-03  1.18e+03  1.20e+03  1.00e+00  9.88e-01  1.86e-06  9.78e-01 33840.3us
%   3:  2.275e-01  2.21e-12  1.15e-06  3.05e-03  1.18e+01  2.65e+01  7.44e-06  4.89e-01  -sigma-   5.38e-01 25190.0us
%   4: -9.293e-01  5.49e-12  5.43e-07  2.02e-03  1.18e-01  1.23e+01  7.44e-06  5.63e-02  -sigma-   2.58e-01 24605.6us
%   5: -3.197e+00  2.29e-10  4.04e-07  3.28e-04  3.37e+00  1.02e+01  7.44e-06  4.32e-01  -sigma-   5.73e-01 24974.5us
%   6: -4.472e+00  1.14e-10  1.77e-07  2.46e-04  3.37e-02  4.48e+00  7.44e-06  1.42e-01  -sigma-   3.57e-01 25210.2us
%   7: -5.377e+00  1.03e-10  1.16e-07  2.01e-04  2.30e-01  2.83e+00  7.44e-06  4.78e-01  -sigma-   7.95e-01 25307.1us
%   8: -6.314e+00  9.29e-11  2.61e-08  9.47e-05  1.52e-02  7.22e-01  7.44e-06  1.34e-01  -sigma-   2.79e-01 24984.0us
%   9: -6.439e+00  7.23e-11  1.92e-08  7.51e-05  8.96e-02  5.09e-01  7.44e-06  6.05e-01  8.32e-02  8.60e-01 25103.3us
%  10: -6.659e+00  3.92e-11  3.43e-09  5.83e-05  3.10e-03  1.10e-01  1.41e-04  2.02e-01  -sigma-   8.86e-01 24927.4us
%  11: -6.699e+00  1.55e-11  5.23e-10  2.40e-05  2.51e-03  4.68e-02  1.41e-04  9.22e-01  1.86e-03  9.80e-01 25145.5us
%  12: -6.737e+00  7.40e-12  1.60e-10  7.29e-07  4.94e-05  1.33e-03  2.89e-07  9.69e-01  4.44e-05  9.91e-01 24892.2us
%  13: -6.738e+00  4.46e-13  6.60e-12  9.80e-09  5.48e-07  1.76e-05  1.67e-08  1.00e+00  5.06e-11  1.00e+00 26835.9us
%  14: -6.738e+00  9.77e-15  5.10e-14  2.92e-09  7.73e-08  5.00e-06  -> clean exit
%  14:status=0x0, cost= -6.73837e+00, |eq|=  5.10e-14, ineq=  2.92e-09,
%               dual=  7.73e-08, gap=  5.00e-06, last alpha=  1.00e+00, |grad|=  9.77e-15 (395494.3us,28249.59us/iter)
% >> [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
% ipmPD_CSsolver.c (skipAffine=0,delta=3,allowSave=0): 300 primal variables, 1 eq. constr., 300 ineq. constr.
% Iter   cost      |grad|      |eq|    inequal     dual      gap       mu      alphaA    sigma     alphaS   time [us]
%  30:<-mx tol->   1.00e-04  1.00e-04                      1.00e-05  1.67e-08
%   1:  1.820e-01  6.13e+02  5.00e-01  1.67e-03  6.00e+02  3.00e+02  1.00e+00  4.94e-01  -sigma-   1.00e+00 28966.2us
%   2:  2.753e-01  4.53e-12  1.20e-06  3.14e-03  1.18e+03  1.20e+03  1.00e+00  9.88e-01  1.86e-06  9.78e-01 25183.4us
%   3:  2.275e-01  2.21e-12  1.15e-06  3.05e-03  1.18e+01  2.65e+01  7.44e-06  4.89e-01  -sigma-   5.38e-01 26771.2us
%   4: -9.293e-01  5.49e-12  5.43e-07  2.02e-03  1.18e-01  1.23e+01  7.44e-06  5.63e-02  -sigma-   2.58e-01 23510.1us
%   5: -3.197e+00  2.29e-10  4.04e-07  3.28e-04  3.37e+00  1.02e+01  7.44e-06  4.32e-01  -sigma-   5.73e-01 26752.8us
%   6: -4.472e+00  1.14e-10  1.77e-07  2.46e-04  3.37e-02  4.48e+00  7.44e-06  1.42e-01  -sigma-   3.57e-01 23427.9us
%   7: -5.377e+00  1.03e-10  1.16e-07  2.01e-04  2.30e-01  2.83e+00  7.44e-06  4.78e-01  -sigma-   7.95e-01 26170.8us
%   8: -6.314e+00  9.29e-11  2.61e-08  9.47e-05  1.52e-02  7.22e-01  7.44e-06  1.34e-01  -sigma-   2.79e-01 23380.8us
%   9: -6.439e+00  7.23e-11  1.92e-08  7.51e-05  8.96e-02  5.09e-01  7.44e-06  6.05e-01  8.32e-02  8.60e-01 26953.6us
%  10: -6.659e+00  3.92e-11  3.43e-09  5.83e-05  3.10e-03  1.10e-01  1.41e-04  2.02e-01  -sigma-   8.86e-01 23647.9us
%  11: -6.699e+00  1.55e-11  5.23e-10  2.40e-05  2.51e-03  4.68e-02  1.41e-04  9.22e-01  1.86e-03  9.80e-01 26729.3us
%  12: -6.737e+00  7.40e-12  1.60e-10  7.29e-07  4.94e-05  1.33e-03  2.89e-07  9.69e-01  4.44e-05  9.91e-01 23330.0us
%  13: -6.738e+00  4.46e-13  6.60e-12  9.80e-09  5.48e-07  1.76e-05  1.67e-08  1.00e+00  5.06e-11  1.00e+00 26453.1us
%  14: -6.738e+00  9.77e-15  5.10e-14  2.92e-09  7.73e-08  5.00e-06  -> clean exit
%  14:status=0x0, cost= -6.73837e+00, |eq|=  5.10e-14, ineq=  2.92e-09,
%               dual=  7.73e-08, gap=  5.00e-06, last alpha=  1.00e+00, |grad|=  9.77e-15 (332535.4us,23752.53us/iter)


% done class2optimizeCS (2.019 sec)
% ipmPD_CSsolver.m (skipAffine=0,delta=3): 300 primal variable, 1 equality constraints, 300 inequality constraints
% Iter   cost      |grad|     |eq|     inequal     dual      gap       mu      alphaA    sigma     alphaS   time [ms]
%  30:<-mx tol->   1.00e-04  1.00e-04                      1.00e-05  1.67e-08
%   1:  1.820e-01  6.13e+02  5.00e-01  1.67e-03  6.00e+02  3.00e+02  1.00e+00  4.94e-01  -sigma-   1.00e+00    71.6ms
%   2:  2.753e-01  9.07e-11  1.20e-06  3.14e-03  1.18e+03  1.20e+03  1.00e+00  9.88e-01  1.86e-06  9.78e-01    60.1ms
%   3:  2.275e-01  3.07e-12  1.15e-06  3.05e-03  1.18e+01  2.65e+01  7.44e-06  4.89e-01  -sigma-   5.38e-01    43.3ms
%   4: -9.293e-01  6.73e-12  5.43e-07  2.02e-03  1.18e-01  1.23e+01  7.44e-06  5.63e-02  -sigma-   2.58e-01    33.8ms
%   5: -3.197e+00  2.30e-10  4.04e-07  3.28e-04  3.37e+00  1.02e+01  7.44e-06  4.32e-01  -sigma-   5.73e-01    57.7ms
%   6: -4.472e+00  1.14e-10  1.77e-07  2.46e-04  3.37e-02  4.48e+00  7.44e-06  1.42e-01  -sigma-   3.57e-01    20.5ms
%   7: -5.377e+00  1.03e-10  1.16e-07  2.01e-04  2.30e-01  2.83e+00  7.44e-06  4.78e-01  -sigma-   7.95e-01    24.9ms
%   8: -6.314e+00  9.29e-11  2.61e-08  9.47e-05  1.52e-02  7.22e-01  7.44e-06  1.34e-01  -sigma-   2.79e-01    30.5ms
%   9: -6.439e+00  7.23e-11  1.92e-08  7.51e-05  8.96e-02  5.09e-01  7.44e-06  6.05e-01  8.32e-02  8.60e-01    36.3ms
%  10: -6.659e+00  3.92e-11  3.43e-09  5.83e-05  3.10e-03  1.10e-01  1.41e-04  2.02e-01  -sigma-   8.86e-01    38.5ms
%  11: -6.699e+00  1.55e-11  5.23e-10  2.40e-05  2.51e-03  4.68e-02  1.41e-04  9.22e-01  1.86e-03  9.80e-01    39.5ms
%  12: -6.737e+00  7.40e-12  1.60e-10  7.29e-07  4.94e-05  1.33e-03  2.89e-07  9.69e-01  4.44e-05  9.91e-01    43.3ms
%  13: -6.738e+00  4.40e-13  6.60e-12  9.80e-09  5.48e-07  1.76e-05  1.67e-08  1.00e+00  5.06e-11  1.00e+00    43.8ms
%  14: -6.738e+00  1.51e-14  5.08e-14  2.92e-09  7.73e-08  5.00e-06  -> clean exit
%  14:status=0x0, cost= -6.73837e+00, |grad|=  1.51e-14, |eq|=  5.08e-14, ineq=  2.92e-09,
%               dual=  7.73e-08, gap=  5.00e-06, last alpha=  1.00e+00 (552.2ms,39.44ms/iter)
% ipmPD_CSsolver.m (skipAffine=0,delta=3): 300 primal variable, 1 equality constraints, 300 inequality constraints
% Iter   cost      |grad|     |eq|     inequal     dual      gap       mu      alphaA    sigma     alphaS   time [ms]
%  30:<-mx tol->   1.00e-04  1.00e-04                      1.00e-05  1.67e-08
%   1:  1.820e-01  6.13e+02  5.00e-01  1.67e-03  6.00e+02  3.00e+02  1.00e+00  4.94e-01  -sigma-   1.00e+00    16.8ms
%   2:  2.753e-01  9.07e-11  1.20e-06  3.14e-03  1.18e+03  1.20e+03  1.00e+00  9.88e-01  1.86e-06  9.78e-01    14.2ms
%   3:  2.275e-01  3.07e-12  1.15e-06  3.05e-03  1.18e+01  2.65e+01  7.44e-06  4.89e-01  -sigma-   5.38e-01    11.2ms
%   4: -9.293e-01  6.73e-12  5.43e-07  2.02e-03  1.18e-01  1.23e+01  7.44e-06  5.63e-02  -sigma-   2.58e-01    15.0ms
%   5: -3.197e+00  2.30e-10  4.04e-07  3.28e-04  3.37e+00  1.02e+01  7.44e-06  4.32e-01  -sigma-   5.73e-01    20.3ms
%   6: -4.472e+00  1.14e-10  1.77e-07  2.46e-04  3.37e-02  4.48e+00  7.44e-06  1.42e-01  -sigma-   3.57e-01    19.1ms
%   7: -5.377e+00  1.03e-10  1.16e-07  2.01e-04  2.30e-01  2.83e+00  7.44e-06  4.78e-01  -sigma-   7.95e-01    26.1ms
%   8: -6.314e+00  9.29e-11  2.61e-08  9.47e-05  1.52e-02  7.22e-01  7.44e-06  1.34e-01  -sigma-   2.79e-01    28.0ms
%   9: -6.439e+00  7.23e-11  1.92e-08  7.51e-05  8.96e-02  5.09e-01  7.44e-06  6.05e-01  8.32e-02  8.60e-01    33.7ms
%  10: -6.659e+00  3.92e-11  3.43e-09  5.83e-05  3.10e-03  1.10e-01  1.41e-04  2.02e-01  -sigma-   8.86e-01    40.9ms
%  11: -6.699e+00  1.55e-11  5.23e-10  2.40e-05  2.51e-03  4.68e-02  1.41e-04  9.22e-01  1.86e-03  9.80e-01    35.0ms
%  12: -6.737e+00  7.40e-12  1.60e-10  7.29e-07  4.94e-05  1.33e-03  2.89e-07  9.69e-01  4.44e-05  9.91e-01    39.0ms
%  13: -6.738e+00  4.40e-13  6.60e-12  9.80e-09  5.48e-07  1.76e-05  1.67e-08  1.00e+00  5.06e-11  1.00e+00    38.8ms
%  14: -6.738e+00  1.51e-14  5.08e-14  2.92e-09  7.73e-08  5.00e-06  -> clean exit
%  14:status=0x0, cost= -6.73837e+00, |grad|=  1.51e-14, |eq|=  5.08e-14, ineq=  2.92e-09,
%               dual=  7.73e-08, gap=  5.00e-06, last alpha=  1.00e+00 (343.1ms,24.51ms/iter)