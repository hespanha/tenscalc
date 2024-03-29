function out=ipmPDminmax_CS(pars)
% See tenscalc/doc/ipm.tex for an explanation of the formulas used here
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha, Raphael Chinchilla).  All rights reserved.

    %packOptimizationVariables=pars.packOptimizationVariables;
    %isSensitivity=pars.isSensitivity;
    addEye2Hessian=pars.addEye2Hessian;

    scaleCost=pars.scaleCost;
    scaleInequalities=pars.scaleInequalities;
    scaleEqualities=pars.scaleEqualities; % currently not used

    cmexfunction=pars.cmexfunction;
    allowSave=pars.allowSave;

    code=pars.code;

    u=pars.u;
    d=pars.d;

    f=pars.objective;

    Fu=pars.Fu;
    Fd=pars.Fd;
    Gu=pars.Gu;
    Gd=pars.Gd;

    nuU=pars.nuU;
    lambdaU=pars.lambdaU;
    nuD=pars.nuD;
    lambdaD=pars.lambdaD;

    nowarningsamesize=true;
    nowarningever=true;

    szHess_=TcheckVariable('Hess_');

    %profile on

    t1=clock();

    %% Define all sizes
    nU=length(u);
    nD=length(d);
    nZ=nU+nD;
    nGu=length(Gu);
    nGd=length(Gd);
    nG=nGu+nGd;
    nFu=length(Fu);
    nFd=length(Fd);
    nF=nFu+nFd;
    fprintf('    # primal vars = %d, # equal constr = %d, # inequal constr = %d...\n',nZ,nG,nF);

    %% Scaling
    if scaleInequalities
        src={};
        dst={};
        if nFu>0
            src{end+1}=abs(1./Fu);
            dst{end+1}=Tvariable('scale4IneqU__',size(Fu));
            Fu=dst{end}.*Fu;
        end
        if nFd>0
            src{end+1}=abs(1./Fd);
            dst{end+1}=Tvariable('scale4IneqD__',size(Fd));
            Fd=dst{end}.*Fd;
        end
        if ~isempty(src)
            declareCopy(code,dst,src,'scaleIneq__');
        end
    end
    if scaleCost>0
        scale4Cost=Tvariable('scale4Cost__',[]);
        declareCopy(code,scale4Cost,abs(scaleCost/f),'scaleCost__');
        f=scale4Cost*f;
        % will also need to scale "desiredGap" to get exist condition that is independent of scaling
        declareGet(code,scale4Cost,'getScale4Cost__');
    end

    fprintf('    getfg()...');
    t2=clock();
    declareGet(code,f,'getf__');

    %% Check minimizer constraints
    Gu_d=gradient(Gu,d);
    if ~strcmp(type(Gu_d),'zeros')
        Gu,Gu_d
        error('cmex2minmaxCS & class2minmaxCS: equality constraints cannot depend on  maximizer optimization variables\n');
    end
    Fu_d=gradient(Fu,d);
    if ~strcmp(type(Fu_d),'zeros')
        Fu,Fu_d
        error('cmex2minmaxCS & class2minmaxCS: equality constraints cannot depend on  maximizer optimization variables\n');
    end

    %% Stack variables
    F=[Fu;Fd];
    lambda=[lambdaU;lambdaD];
    G=Tzeros(0);
    G=[Gu;Gd];
    nu=[nuU;nuD];
    gap=lambda*F;

    %% Declare gets for exit condition and output
    if nF>0
        mu=Tvariable('mu__',[],nowarningsamesize,nowarningever);

        declareSet(code,mu,'setMu__');
        declareGet(code,{gap,min(F,[],1),min(lambda,[],1)},'getGapMinFMinLambda__');
    end
    if nG>0
        declareGet(code,norminf(G),'getNorminf_G__');
    end

    %% Automatic initialization of the lambda variable
    dst={};
    src={};
    if nGu>0
        dst{end+1}=nuU;
        src{end+1}=Tones(nGu);
    end
    if nGd>0
        dst{end+1}=nuD;
        src{end+1}=Tones(nGd);
    end
    if ~isempty(dst)
        declareCopy(code,dst,src,'initDualEq__');
    end
    if nF>0
        dst={};
        src={};
        if nFu>0
            dst{end+1}=lambdaU;
            src{end+1}=mu*Tones(nFu)./Fu;
        end
        if nFd>0
            dst{end+1}=lambdaD;
            src{end+1}=mu*Tones(nFd)./Fd;
        end
        declareCopy(code,dst,src,'initDualIneq__');
    end

    %% Construct the "Hessian" matrix
    Lf=f;
    if nGu>0
        Lf=Lf+nuU*Gu;
    end
    if nGd>0
        Lf=Lf+nuD*Gd;
    end
    if nFu>0
        Lf=Lf-lambdaU*Fu;
    end
    if nFd>0
        Lf=Lf+lambdaD*Fd;
    end

    alphaPrimal=Tvariable('alphaPrimal__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaPrimal,'setAlphaPrimal__');
    alphaDualEq=Tvariable('alphaDualEq__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaDualEq,'setAlphaDualEq__');
    alphaDualIneq=Tvariable('alphaDualIneq__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaDualIneq,'setAlphaDualIneq__');

    fprintf('(%.2f sec)\n    1st derivates...',etime(clock(),t2));
    t2=clock();
    if addEye2Hessian
        addEye2HessianU=Tvariable('addEye2HessianU__',[],nowarningsamesize,nowarningever);
        addEye2HessianD=Tvariable('addEye2HessianD__',[],nowarningsamesize,nowarningever);
        addEye2HessianEq=Tvariable('addEye2HessianEq__',[],nowarningsamesize,nowarningever);
        declareSet(code,addEye2HessianU,'setAddEye2HessianU__');
        declareSet(code,addEye2HessianD,'setAddEye2HessianD__');
        declareSet(code,addEye2HessianEq,'setAddEye2HessianEq__');
    else
        addEye2HessianU=Tzeros([]);
        addEye2HessianD=Tzeros([]);
        addEye2HessianEq=Tzeros([]);
    end
    out.addEye2HessianU=addEye2HessianU;
    out.addEye2HessianD=addEye2HessianD;
    out.addEye2HessianEq=addEye2HessianEq;

    Lf_z=gradientVector(Lf,{u,d});
    % for exit condition
    declareGet(code,norminf(Lf_z),'getNorminf_Grad__');

    fprintf('(%.2f sec)\n    2nd derivatives...',etime(clock(),t2));
    t2=clock();

    % Derivatives needed to compute WW 
    Lf_zz=gradientVector(Lf_z,{u,d});
    Lf_zl=gradientVector(Lf_z,{lambdaU,lambdaD});
    Lf_zn=gradientVector(Lf_z,{nuU,nuD});

    % Derivatives needed to compute inertial
    Lf_u=gradient(Lf,u);
    Lf_d=gradient(Lf,d);
    Lf_ud=gradient(Lf_u,d);
    Lf_uu=gradient(Lf_u,u);
    Lf_dd=gradient(Lf_d,d);
    Lf_du=gradient(Lf_d,u);

    fprintf('(%.2f sec)\n    WW & b...',etime(clock(),t2));
    t2=clock();

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Hessian computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    WWD = Lf_dd-addEye2HessianD*Teye([nD,nD]);
    WWUD=[Lf_uu+addEye2HessianU*Teye([nU,nU]),Lf_ud;
          Lf_du,                              WWD];
    if nF>0
        WW=[WWUD,Lf_zn,Lf_zl;
            gradientVector([Gu,Gd],{u,d}),-addEye2HessianEq*Teye([nG,nG]),Tzeros([nG,nF]);
            gradientVector([-Fu;Fd],{u,d}),Tzeros([nF,nG]),diag([-Fu./lambdaU;Fd./lambdaD])];
    else
        WW=[WWUD,Lf_zn;
            gradientVector([Gu,Gd],{u,d}),-addEye2HessianEq*Teye([nG,nG])];
    end
    out.Lf_z=Lf_z;
    out.Lf_zz=Lf_zz;
    out.Hess=WW;

    if nF>0
        b_s=[-Lf_z;
             -G;
             +Fu-mu./lambdaU;
             -Fd+mu./lambdaD];
    else
        b_s=[-Lf_z;
             -G];
    end

    WW=declareAlias(code,WW,'WW__',false,nowarningsamesize,nowarningever);

    ldl_ww=ldl(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Inertial computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    tol=0*1e-8; % minimum eigenvalue to be considered positive
    
    if nF>0
        HessD=[Lf_dd-addEye2HessianD*Teye([nD,nD]),gradientVector(Lf_d,{nuD}),gradientVector(Lf_d,{lambdaD});
               gradient(Gd,d),-addEye2HessianEq*Teye([nGd,nGd]),Tzeros([nGd,nFd]);
               gradient(Fd,d),Tzeros([nFd,nGd]),diag(Fd./lambdaD)];
    else
        HessD=[Lf_dd-addEye2HessianD*Teye([nD,nD]),gradientVector(Lf_d,{nuD});
               gradient(Gd,d),-addEye2HessianEq*Teye([nGd,nGd])];
    end
    out.HessD=HessD;    
    out.d_HessD=ldl_d(ldl(HessD));
    declareGet(code,{sum(heaviside(out.d_HessD-tol)),...
                     sum(heaviside(-out.d_HessD-tol))},'getHessDinertia__');

    out.HessU=WW;
    out.d_HessU=ldl_d(ldl_ww);
    declareGet(code,{sum(heaviside(out.d_HessU-tol)),...
                     sum(heaviside(-out.d_HessU-tol))},'getHessUinertia__');
    
    fprintf('(%.2f sec)\n    adding WW...',etime(clock(),t2));
    t2=clock();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% search direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('(%.2f sec)\n    combined direction...',etime(clock(),t2));
    t2=clock();

    dx_s=ldl_ww\b_s;
    dx_s=declareAlias(code,dx_s,'dx_s__',false,nowarningsamesize,nowarningever);

    declareGet(code,{norminf((WW*dx_s-b_s))},'getDirectionError__');

    dU_s=dx_s(1:nU);
    dD_s=dx_s(nU+1:nU+nD);
    dZ_s=dx_s(1:nZ);
    dNu_s=dx_s(nZ+1:nZ+nG);
    dNuU_s=dx_s(nZ+1:nZ+nGu);
    dNuD_s=dx_s(nZ+nGu+1:nZ+nGu+nGd);

    if nF>0
        dLambda_s=dx_s(nZ+nG+1:nZ+nG+nF);
    else
        dLambda_s=Tzeros(nF);
    end

    fprintf('(%.2f sec)\n    new vectores...',etime(clock(),t2));
    t2=clock();

    newU_s=u+alphaPrimal*dU_s;
    newD_s=d+alphaPrimal*dD_s;

    if nF>0
        F_z=gradientVector(F,{u,d});
        maxAlphaPrimal_s=clp(F,F_z*dZ_s);
        maxAlphaDualIneq_s=clp(lambda,dLambda_s);
        newLambda_s=lambda+alphaDualIneq*dLambda_s;
        newF_s=substitute(F,u,newU_s);
        newF_s=substitute(newF_s,d,newD_s);
        declareGet(code,{maxAlphaPrimal_s,maxAlphaDualIneq_s},'getMaxAlphas_s__');
        declareGet(code,min(newF_s,[],1),'getMinF_s__');
    end
    
    fprintf('(%.2f sec)\n    updatePrimalDual()...',etime(clock(),t2));
    t2=clock();
    dst={u,d};
    src={newU_s,newD_s};
    if nGu>0
        dst{end+1}=nuU;
        src{end+1}=nuU+alphaDualEq*dNuU_s;
    end
    if nGd>0
        dst{end+1}=nuD;
        src{end+1}=nuD+alphaDualEq*dNuD_s;
    end
    if nFu>0
        dst{end+1}=lambdaU;
        src{end+1}=newLambda_s(1:nFu);
    end
    if nFd>0
        dst{end+1}=lambdaD;
        src{end+1}=newLambda_s(nFu+1:nF);
    end
    
    declareCopy(code,dst,src,'updatePrimalDual__');
    
    if ~isempty(szHess_) && ~myisequal(szHess_,size(out.Hess))
        warning('\nvariable: ''Hess_'' already exists with the wrong size [%d,%d], should be [%d,%d]\n',...
                szHess_(1),szHess_(2),size(out.Hess,1),size(out.Hess,2));
    end
    
    fprintf('(%.2f sec)\n    ',etime(clock(),t2));
    fprintf(' done ipmPDminmax_CS symbolic computations (%.3f sec)\n',etime(clock(),t1));

    %profile off
    %profile viewer

end
