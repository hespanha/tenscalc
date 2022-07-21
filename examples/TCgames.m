% This script uses Tenscals's function `cmex2equilibriumLatentCS` to
% generate solvers for several zero-sum sum games.
%
% Solvers can be generated once and then called multiple
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

help TCgames

createSolvers=input('Do you want to create the solvers (Y/n)? ','s');
createSolvers=~isequal(createSolvers,'n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero-sum matrix-game
%% Computation of mixed saddle-point
%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%%%%
if true
    N1=50;
    N2=20;

    if createSolvers

        %% Generate solver
        Tcalculus.clear();

        Tvariable A1 [N1,N2];
        Tvariable u N1;   % mixed policy for P1
        Tvariable d N2;   % mixed policy for P2

        J1=u*A1*d;   % transpose on 1st u not needed: u is a vactor => u*A1 is also a vector => (u*A1)*d is the inner product of two vectors

        fprintf('%dx%d matrix game saddle point\n',N1,N2);
        classname=cmex2equilibriumLatentCS('classname','tmp_matrixsaddle',...
                                           'P1objective',J1,...
                                           'P2objective',-J1,...
                                           'P1optimizationVariables',{u},...
                                           'P2optimizationVariables',{d},...
                                           'P1constraints',{sum(u,1)==1,u>=0},...  % mixed policy constraint
                                           'P2constraints',{sum(d,1)==1,d>=0},...  % mixed policy constraint
                                           'outputExpressions',{u,d,J1},...
                                           'parameters',{A1},...
                                           'compilerOptimization','-O0',...
                                           'solverVerboseLevel',2);
    end

    %% Use solver
    obj=tmp_matrixsaddle();

    s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
    RandStream.setGlobalStream(s);

    % Set optimization parameters
    thisA1=rand(N1,N2);
    setP_A1(obj,thisA1);

    % Initialize primal variables
    u0=ones(N1,1)/N1;
    d0=ones(N2,1)/N2;
    setV_u(obj,u0);
    setV_d(obj,d0);
    % Solve optimization
    mu0=1;
    maxIter=30;
    saveIter=-1;

    % Call solver
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [u1,d1,j1]=getOutputs(obj);

    j1

    if false
        disp('Verification using linprog')
        [v2,u2,d2]=saddle_mnmx(thisA1);
        tic
        [v2,u2,d2]=saddle_mnmx(thisA1);
        toc
        v2
        %u2'
        %d2'
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero-sum quadratic cost
%% Computation of saddle-point equlibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    N=200;

    if createSolvers

        %% Generate solver
        Tcalculus.clear();

        Tvariable u N;  % policy for P1
        Tvariable d N;  % policy for P2

        Ju=norm2(u+d-2)-10*norm2(d);

        fprintf('(%d,%d) quadratic saddle-point\n',N,N);
        classname=cmex2equilibriumLatentCS('classname','tmp_qsaddle',...
                                           'P1objective',Ju,...
                                           'P2objective',-Ju,...
                                           'P1optimizationVariables',{u},...
                                           'P2optimizationVariables',{d},...
                                           'P1constraints',{u>=-1,u<=1},...  % bound on absolute value
                                           'P2constraints',{d>=-1,d<=1},...  % bound on absolute value
                                           'outputExpressions',{u,d,Ju},...
                                           'compilerOptimization','-O0',...
                                           'solverVerboseLevel',2);
    end

    %% Use solver
    obj=tmp_qsaddle();

    s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
    RandStream.setGlobalStream(s);

    % Initialize primal variables
    u0=.5*ones(N,1);
    d0=.5*ones(N,1);
    setV_u(obj,u0);
    setV_d(obj,d0);
    % Solve optimization
    mu0=1;
    maxIter=30;
    saveIter=-1;
    % Call solver
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [u1,d1,ju]=getOutputs(obj);

    ju
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero-sum quadratic cost, with latent (state) variable
%% Computation of saddle-point equlibrium
%%
%% Note: This problem is exactly the same as the previous
%%       one, but formulated using a latent variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    N=200;

    if createSolvers

        %% Generate solver
        Tcalculus.clear();

        Tvariable u N;   % policy for P1
        Tvariable d N;   % policy for P2
        Tvariable x N;   % latent variables (will be constrained to be u+d)

        Ju=norm2(x-2)-10*norm2(d);
        Jd=norm2(x-2)-10*norm2(d);

        fprintf('(%d,%d) quadratic saddle-point with state\n',N,N);
        classname=cmex2equilibriumLatentCS('classname','tmp_qsaddlestate',...
                                           'P1objective',Ju,...
                                           'P2objective',-Jd,...
                                           'P1optimizationVariables',{u},...
                                           'P2optimizationVariables',{d},...
                                           'latentVariables',{x},...
                                           'P1constraints',{u>=-1,u<=1},...  % bound on absolute value
                                           'P2constraints',{d>=-1,d<=1},...  % bound on absolute value
                                           'latentConstraints',{x==u+d},...  % constraint on latent variable
                                           'outputExpressions',{u,d,x,Ju,Jd},...
                                           'compilerOptimization','-O0',...
                                           'solverVerboseLevel',2);
    end

    %% Use solver
    s = RandStream('mt19937ar','Seed',0); % fix seed for reproducibility
    RandStream.setGlobalStream(s);

    obj=tmp_qsaddlestate();

    % Initialize primal variables
    u0=zeros(N,1);
    d0=zeros(N,1);
    x0=u0+d0;
    setV_u(obj,u0);
    setV_d(obj,d0);
    setV_x(obj,x0);
    % Solve optimization
    mu0=1;
    maxIter=30;
    saveIter=-1;
    % Call solver
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [u1,d1,x1,ju,jd]=getOutputs(obj);


    ju,jd
end