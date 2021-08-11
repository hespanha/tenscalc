% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

% Many ways to skin a cat...
%
% This example solve several 2-player games.

clear all;
% remove previous solvers
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Matrix-game Saddle-point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%
if 1
    clear all;

    N1=50;
    N2=20;

    Tvariable A1 [N1,N2];
    Tvariable u N1;
    Tvariable d N2;

    J1=tprod(u,[-1],A1,[-1,-2],d,[-2]);

    fprintf('%dx%d matrix game saddle point\n',N1,N2);
    classname=cmex2equilibriumLatentCS('classname','tmp_matrixsaddle',...
                                       'P1objective',J1,...
                                       'P2objective',-J1,...
                                       'P1optimizationVariables',{u},...
                                       'P2optimizationVariables',{d},...
                                       'P1constraints',{sum(u,1)==1,u>=0},...
                                       'P2constraints',{sum(d,1)==1,d>=0},...
                                       'outputExpressions',{u,d,J1},...
                                       'parameters',{A1},...
                                       'compilerOptimization','-O0',...
                                       'solverVerboseLevel',3);

    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);

    obj=feval(classname);

    thisA1=rand(N1,N2);

    % Set parameters
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
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [u1,d1,j1]=getOutputs(obj);

    j1
    %u1'
    %d1'
    status

    if 0
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 2-Matrix game Nash
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    clear all;

    N1=50;
    N2=20;

    Tvariable A1 [N1,N2];
    Tvariable A2 [N1,N2];
    Tvariable u N1;
    Tvariable d N2;

    J1=tprod(u,[-1],A1,[-1,-2],d,[-2]);
    J2=tprod(u,[-1],A2,[-1,-2],d,[-2]);

    fprintf('%dx%d matrix game nash equilibrium\n',N1,N2);
    classname=cmex2equilibriumLatentCS('classname','tmp_matrixnash',...
                                       'P1objective',J1,...
                                       'P2objective',J2,...
                                       'P1optimizationVariables',{u},...
                                       'P2optimizationVariables',{d},...
                                       'P1constraints',{sum(u,1)==1,u>=0},...
                                       'P2constraints',{sum(d,1)==1,d>=0},...
                                       'outputExpressions',{u,d,J1,J2},...
                                       'parameters',{A1,A2},...
                                       'compilerOptimization','-O0',...
                                       'solverVerboseLevel',3);

    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);

    obj=feval(classname);

    thisA1=rand(N1,N2);
    thisA2=-thisA1;
    thisA2=rand(N1,N2);

    % Set parameters
    setP_A1(obj,thisA1);
    setP_A2(obj,thisA2);
    % Initialize primal variables
    u0=ones(N1,1)/N1;
    d0=ones(N2,1)/N2;
    setV_u(obj,u0);
    setV_d(obj,d0);
    % Solve optimization
    mu0=1;
    maxIter=30;
    saveIter=-1;
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [u1,d1,j1,j2]=getOutputs(obj);

    j1,j2
    %u1'
    %d1'
    status

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Quadratic saddle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    clear all;
    N=200

    Tvariable u N;
    Tvariable d N;

    Ju=norm2(u+d-2)-10*norm2(d);

    fprintf('(%d,%d) quadratic saddle-point\n',N,N);
    classname=cmex2equilibriumLatentCS('classname','tmp_qsaddle',...
                                       'P1objective',Ju,...
                                       'P2objective',-Ju,...
                                       'P1optimizationVariables',{u},...
                                       'P2optimizationVariables',{d},...
                                       'P1constraints',{u>=-1,u<=1},...
                                       'P2constraints',{d>=-1,d<=1},...
                                       'outputExpressions',{u,d,Ju},...
                                       'compilerOptimization','-O0',...
                                       'solverVerboseLevel',3);

    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);

    obj=feval(classname);

    % Initialize primal variables
    u0=.5*ones(N,1);
    d0=.5*ones(N,1);
    setV_u(obj,u0);
    setV_d(obj,d0);
    % Solve optimization
    mu0=1;
    maxIter=30;
    saveIter=-1;
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [u1,d1,ju]=getOutputs(obj);

    ju
    status
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Quadratic saddle with state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    clear all;
    N=200

    Tvariable u N;
    Tvariable d N;
    Tvariable xu N;
    Tvariable xd N;

    Ju=norm2(xu-2)-10*norm2(d);
    Jd=norm2(xd-2)-10*norm2(d);

    fprintf('(%d,%d) quadratic saddle-point with state\n',N,N);
    classname=cmex2equilibriumLatentCS('classname','tmp_qsaddlestate',...
                                       'P1objective',Ju,...
                                       'P2objective',-Jd,...
                                       'P1optimizationVariables',{u,xu},...
                                       'P2optimizationVariables',{d,xd},...
                                       'P1constraints',{xu==u+d,u>=-1,u<=1},...
                                       'P2constraints',{xd==u+d,d>=-1,d<=1},...
                                       'outputExpressions',{u,xu,d,xd,Ju,Jd},...
                                       'compilerOptimization','-O0',...
                                       'solverVerboseLevel',3);
    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);

    obj=feval(classname);

    % Initialize primal variables
    u0=zeros(N,1);
    d0=zeros(N,1);
    xu0=u0+d0;
    xd0=u0+d0;
    setV_u(obj,u0);
    setV_xu(obj,xu0);
    setV_d(obj,d0);
    setV_xd(obj,xd0);
    % Solve optimization
    mu0=1;
    maxIter=30;
    saveIter=-1;
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    [status,iter,time]=solve(obj,mu0,int32(maxIter),int32(saveIter));
    % Get outputs
    [u1,xu1,d1,xd1,ju,jd]=getOutputs(obj);
    [u1,xu1,d1,xd1,ju,jd]=getOutputs(obj);

    ju,jd
    norm(xu1-xd1)
    status
end