function out=ipmPD_CStimesLambda(pars)
% See tenscalc/doc/ipm.tex for an explanation of the formulas used here
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

%% Hardcoded parameters
nowarningsamesize=true;
nowarningever=true;

%% Retrieve input paramneters
packOptimizationVariables=pars.packOptimizationVariables;
isSensitivity=pars.isSensitivity;
smallerNewtonMatrix=pars.smallerNewtonMatrix;
skipAffine=pars.skipAffine;

addEye2Hessian=pars.addEye2Hessian;
scaleCost=pars.scaleCost;
scaleInequalities=pars.scaleInequalities;
scaleEqualities=pars.scaleEqualities; % currently not used

useLDL=pars.useLDL;
atomicFactorization=pars.atomicFactorization;

cmexfunction=pars.cmexfunction;
allowSave=pars.allowSave;
debugConvergence=pars.debugConvergence;

if any(isSensitivity)
    error('ipmPD_CStimesLambda: isSensitivity=true not implemented, use ipmPD_CS.m \n');
end
if smallerNewtonMatrix==true
    error('ipmPD_CStimesLambda: smallerNewtonMatrix=true not implemented, use ipmPD_CS.m \n');
end
if skipAffine==false
    error('ipmPD_CStimesLambda: skipAffine=false not implemented, use ipmPD_CS.m \n');
end
if useLDL==false
    error('ipmPD_CStimesLambda: useLDL=false not implemented, use ipmPD_CS.m \n');
end
if atomicFactorization==true
    error('ipmPD_CStimesLambda: atomicFactorization=true not implemented, use ipmPD_CS.m \n');
end
if debugConvergence==true
    error('ipmPD_CStimesLambda: debugConvergence=true not implemented, use ipmPD_CS.m \n');
end

code=pars.code;

%% Retrieve primal variables, dual variables, and constraints
u=pars.u;
f=pars.f;
F=pars.F;
G=pars.G;

if ~packOptimizationVariables
    uList=pars.uList;
    f_uList=pars.f_uList;
    F_uList=pars.F_uList;
    G_uList=pars.G_uList;
end

lambda=pars.lambda;
nu=pars.nu;

out.u=u;
out.F=F;
out.G=G;
out.nu=nu;
out.lambda=lambda;

%% Define all sizes
szHess_=TcheckVariable('Hess_');

nU=length(u);
nG=length(G);
nF=length(F);
fprintf('    # primal vars = %d, # equal constr = %d, # inequal constr = %d...\n',nU,nG,nF);

fprintf('\n  Starting ipmPD_CStimesLambda symbolic computations...\n');
t1=clock();

%% Scaling
if nF>0 && scaleInequalities
    scale4Ineq=Tvariable('scale4Ineq__',size(F));
    declareCopy(code,scale4Ineq,abs(1./F),'scaleIneq__');
    F=scale4Ineq.*F;
    if ~packOptimizationVariables
        F_uList=scale4Ineq.*F_uList;
    end
end
if scaleCost>0
    scale4Cost=Tvariable('scale4Cost__',[]);
    declareCopy(code,scale4Cost,abs(scaleCost/f),'scaleCost__');
    f=scale4Cost*f;
    if ~packOptimizationVariables
        f_uList=scale4Cost*f_uList;
    end
    % will also need to scale "desiredGap" to get exist condition that is independent of scaling
    declareGet(code,scale4Cost,'getScale4Cost__');
end

%% Matrices used to modify Hessian's inertia
if addEye2Hessian
    addEye2HessianU=Tvariable('addEye2HessianU__',[],nowarningsamesize,nowarningever);
    addEye2HessianEq=Tvariable('addEye2HessianEq__',[],nowarningsamesize,nowarningever);
    declareSet(code,addEye2HessianU,'setAddEye2HessianU__');
    declareSet(code,addEye2HessianEq,'setAddEye2HessianEq__');
else
    addEye2HessianU=Tzeros([]);
    addEye2HessianEq=Tzeros([]);
end
out.addEye2HessianU=addEye2HessianU;
out.addEye2HessianEq=addEye2HessianEq;

%% Cost
fprintf('    getJ()...');
t2=clock();
declareGet(code,f,'getJ__');


%% Compute Lagrangian and 1st derivatives
fprintf('(%.2f sec)\n    1st derivatives...',etime(clock(),t2));
t2=clock();
if packOptimizationVariables
    f_u=gradient(f,u);
    Lf=f;
else
    f_u=gradientVector(f_uList,uList);
    Lf=f_uList;
end
out.Lf_u=f_u;
% inequality constraints
if nF>0
    if packOptimizationVariables
        F_u=gradient(F,u);
        gap=tprod(lambda,-1,F,-1);                        % gap=lambda*F;
    else
        F_u=gradientVector(F_uList,uList);
        gap=tprod(lambda,-1,F_uList,-1);                  % gap=lambda*F;
    end
    Lf=Lf-gap;                                            % Lf=Lf-gap;
    % out.Lf_u=out.Lf_u-tprod(F_u,[-1,1],lambda,-1);        % Lf_u=Lf_u-F_u'*lambda;
else
    F=Tzeros(0);
    F_u=Tzeros([0,nU]);
    gap=Tzeros([]);
end
% equality constraints
if nG>0
    if packOptimizationVariables
        G_u=gradient(G,u);
        Lf=Lf+tprod(nu,-1,G,-1);                              % Lf=Lf+nu*G;
    else
        G_u=gradientVector(G_uList,uList);
        Lf=Lf+tprod(nu,-1,G_uList,-1);                        % Lf=Lf+nu*G;
    end
    % out.Lf_u=out.Lf_u+tprod(G_u,[-1,1],nu,-1);                % Lf_u=Lf_u+G_u'*nu;
else
    G=Tzeros(0);
    G_u=Tzeros([0,nU]);
end

%% Compute Lagrangian's 1st and 2nd derivative
fprintf('(%.2f sec)\n    2nd derivatives...',etime(clock(),t2));
t2=clock();
if packOptimizationVariables
    %out.Lf_uu=gradient(out.Lf_u,u);
    [out.Lf_uu,out.Lf_u]=hessian(Lf,u); % seems to be more efficient than using the previously computed F_u/G_u
else
    [out.Lf_u,out.Lf_uu]=gradientVector(Lf,uList); % seems to be more efficient than using the previously computed F_u/G_u

    % Replace uList by u
    Lf=substitute(Lf,uList,u);
    f_u=substitute(f_u,uList,u);
    out.Lf_u=substitute(out.Lf_u,uList,u);
    F_u=substitute(F_u,uList,u);
    gap=substitute(gap,uList,u);
    G_u=substitute(G_u,uList,u);
    out.Lf_uu=substitute(out.Lf_uu,uList,u);
end

% gets for Lagrangian's derivatives
declareGet(code,out.Lf_u,'getGrad__');
declareGet(code,norminf(out.Lf_u),'getNorminf_Grad__');

%% Barrier variable
if nF>0
    out.mu=Tvariable('mu__',[],nowarningsamesize,nowarningever);
    muOnes=reshape(out.mu,1);
    muOnes=repmat(muOnes,nF);

    declareSet(code,out.mu,'setMu__');

    % Automatic initialization of lambda
    declareCopy(code,lambda,muOnes./F,'initDualIneq__');

    %% Current value for primal and dual inequalities and gap
    declareGet(code,{gap,min(F,[],1),min(lambda,[],1)},'getGapMinFMinLambda__');
else
    out.mu=Tzeros([]);
    muOnes=Tzeros(0);
end

if nG>0
    % Automatic initialization of nu
    % simple: just set to 1
    declareCopy(code,nu,Tones(nG),'initDualEq__');
    % complex: least squares
    WW0=[Teye(nU,nU),G_u';G_u,-addEye2HessianEq*Teye(nG,nG)];
    factor_ww0=ldl(WW0,[cmexfunction,'_WW0.subscripts'],[cmexfunction,'_WW0.values']);
    b0=[F_u'*lambda-f_u;Tzeros(nG)];
    wnu0=factor_ww0\b0;
    wnu0=full(wnu0(nU+1:end)); % can easily get sparse when cost is linear on opt. variables
    declareCopy(code,nu,wnu0(nU+1:end),'initDualEqX__');

    %% Error in equality constraints
    declareGet(code,norminf(G),'getNorminf_G__');
end

%% Step size
alphaPrimal=Tvariable('alphaPrimal__',[],nowarningsamesize,nowarningever);
declareSet(code,alphaPrimal,'setAlphaPrimal__');
alphaDualEq=Tvariable('alphaDualEq__',[],nowarningsamesize,nowarningever);
declareSet(code,alphaDualEq,'setAlphaDualEq__');
alphaDualIneq=Tvariable('alphaDualIneq__',[],nowarningsamesize,nowarningever);
declareSet(code,alphaDualIneq,'setAlphaDualIneq__');

fprintf('(%.2f sec)\n    WW...',etime(clock(),t2));
t2=clock();

%% Compute Newton step
WW11=out.Lf_uu+addEye2HessianU*Teye(size(out.Lf_uu));
WW=[WW11,G_u',-F_u'*diag(lambda);
    G_u,-addEye2HessianEq*Teye([nG,nG]),Tzeros([nG,nF]);
    -diag(lambda)*F_u,Tzeros([nF,nG]),-diag(F.*lambda)];

out.Hess=WW;

factor_ww=ldl(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);

out.dHess=ldl_d(factor_ww);
out.lHess=ldl_l(factor_ww);
tol=0*1e-8; % minimum eigenvalue to be considered positive -- as we let addEye2Hessian get smaller, this generally also needs to get smaller
declareGet(code,{sum(heaviside(out.dHess-tol)),...
    sum(heaviside(-out.dHess-tol))},'getHessInertia__');

b=[-out.Lf_u;
    -G;
    lambda.*F-muOnes];

%% Newton direction
dx=factor_ww\b;
dx=declareAlias(code,dx,'dx__',false,nowarningsamesize,nowarningever);

declareGet(code,{norminf((WW*dx-b))},'getDirectionError__');

%% Update primal and dual variables
dU=dx(1:nU);
newU=u+alphaPrimal*dU;
dNu=dx(nU+1:nU+nG);
newNu=nu+alphaDualEq*dNu;

curvature=dU*(WW11*dU);
declareGet(code,curvature,'getCurvature__');

if nF>0
    dLambda=dx(nU+nG+1:nU+nG+nF);
    newLambda=lambda.*(Tones(nF)+alphaDualIneq*dLambda);

    maxAlphaPrimal=clp(F,F_u*dU);
    maxAlphaDualIneq=clp(Tones(nF),dLambda);
    declareGet(code,{maxAlphaPrimal,maxAlphaDualIneq},'getMaxAlphas_s__');
    newF=substitute(F,u,newU);
    declareGet(code,min(newF,[],1),'getMinF_s__');
else
    newLambda=Tzeros(nF);
end
declareCopy(code,{u,nu,lambda},{newU,newNu,newLambda},'updatePrimalDual__');

out.b=b;
out.dx=dx;

% declareSave after lu's to make sure previous "typical" values
% are used by lu, prior to being overwritten by declareSave
if allowSave
    declareSave(code,WW,'saveWW__',[cmexfunction,'_WW.subscripts'])
end

if ~isempty(szHess_) && ~myisequal(szHess_,size(out.Hess))
    warning('\nvariable: ''Hess_'' already exists with the wrong size [%d,%d], should be [%d,%d]\n',...
        szHess_(1),szHess_(2),size(out.Hess,1),size(out.Hess,2));
end

fprintf('(%.2f sec)\n    ',etime(clock(),t2));
fprintf('  done ipmPD_CStimesLambda symbolic computations (%.3f sec)\n',etime(clock(),t1));

%    profile off;profile viewer

end
