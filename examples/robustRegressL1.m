% This script uses Tenscals's function `cmex2equilibriumLatentCS` to
% generate a solver for an L1-regularized regression of the form:
%
%    minimize || y - theta0 ones(size(y)) - H \theta ||_2 + lambda ||\theta||_1
%
%    w.r.t. theta0 in R, theta in R^n
%
% The solver can be generated once and then called multiple
% times. Creating a solver typically takes much longer than solving
% the optimization by creating the solver.
%
% This script give you the option and skipping creating the solver,
% but you can only skip generating a solver if you have done it once.

% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all;
% remove previous solvers
%delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

help robustRegressL1

createSolvers=input('Do you want to create the solvers (Y/n)? ','s');
createSolvers=~isequal(createSolvers,'n');

%%%%%%%%%%%%%%%%%%%
%% Load test data
%%%%%%%%%%%%%%%%%%%

m=1000; % size of training data
n=15;  % # of independent variables

if createSolvers

    %%%%%%%%%%%%%%%%%%%
    %% Generate solver
    %%%%%%%%%%%%%%%%%%%

    Tcalculus.clear();

    Tvariable lambda [];

    Tvariable theta0   [];
    Tvariable theta    [n];
    Tvariable absTheta [n];
    Tvariable v    [];

    if 1
        Tvariable y [m];
        Tvariable H [m,n];

        v2=norm2(y-theta0*Tones(m)-H*theta);
    else
        Tvariable y2   [];
        Tvariable sumY [];
        Tvariable yH   [n];
        Tvariable sumH [n];
        Tvariable HH   [n,n];
        v2=y2-2*sumY*theta0-2*yH*theta+2*theta0*sumH*theta+theta0*theta0*m+theta*HH*theta;
    end

    v=sqrt(v2);
    J =v+lambda*sum(absTheta,1);

    classname=cmex2optimizeCS('classname','tmp_lasso',...
                              'objective',J,...
                              'optimizationVariables',{theta0,theta,absTheta},...
                              'constraints',{absTheta>theta;absTheta>-theta},...
                              'parameters',{lambda,y,H},...
                              'outputExpressions',{theta,theta0,J},...
                              'scaleCost',1,...
                              'solverVerboseLevel',2);


    preMultiplied=true;
    useSqrt=true;      % true is slightly slower, but can be more robust
    smoothSqrt=false;  % true generally leads to slower convergence

end

%%%%%%%%%%%%%%%
%% Use solver
%%%%%%%%%%%%%%%

% Create object
obj=tmp_lasso();

s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
RandStream.setGlobalStream(s);

lambda=10;
setP_lambda(obj,lambda);

for i=1:3

    %% Generate synthetic training data
    thisTheta=randn(n,1);
    thisTheta(rand(n,1)<.5)=0;
    thisTheta0=randn(1,1);
    thisH=randn(m,n);
    thisy=thisTheta0+thisH*thisTheta+.2*randn(m,1);

    %% Set parameters
    setP_y(obj,thisy);
    setP_H(obj,thisH);

    %% Random initial value for solver
    theta=.01*randn(n,1);
    theta0=mean(thisy-thisH*theta);

    setV_theta0(obj,theta0);
    setV_theta(obj,theta);
    setV_absTheta(obj,abs(theta)+1);

    % Solve optimization
    mu0=.001;
    maxIter=100;
    saveIter=-1;
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [theta,theta0,J]=getOutputs(obj);

    fprintf('Cost = %.2f, theta0 = %.3f, theta =\n',J,theta0);
    disp(theta');

    if status
        error()
    end
end
