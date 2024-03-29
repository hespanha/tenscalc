function out=ipmPD_CS(pars)
% See tenscalc/doc/ipm.tex for an explanation of the formulas used here
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    %% Retrieve input paramneters
    packOptimizationVariables=pars.packOptimizationVariables;
    isSensitivity=pars.isSensitivity;
    smallerNewtonMatrix=pars.smallerNewtonMatrix;
    addEye2Hessian=pars.addEye2Hessian;
    skipAffine=pars.skipAffine;

    scaleCost=pars.scaleCost;
    scaleInequalities=pars.scaleInequalities;
    scaleEqualities=pars.scaleEqualities; % currently not used

    useLDL=pars.useLDL;
    atomicFactorization=pars.atomicFactorization;

    cmexfunction=pars.cmexfunction;
    allowSave=pars.allowSave;
    debugConvergence=pars.debugConvergence;

    code=pars.code;

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


    %profile clear;profile on

    nowarningsamesize=true;
    nowarningever=true;

    trustRegion=false;

    szHess_=TcheckVariable('Hess_');

    if smallerNewtonMatrix
        fprintf('\n  Starting ipmPD_CS symbolic computations (smallNewtonMatrix)...\n');
    else
        fprintf('\n  Starting ipmPD_CS symbolic computations (largeNewtonMatrix)...\n');
    end
    t1=clock();

    %% Select factorization method
    if useLDL
        if atomicFactorization
            factor=@lu_sym;
            error('Bug, should not come here');
        else
            factor=@ldl;
        end
    else
        if atomicFactorization
            factor=@lu_sym;
        else
            factor=@lu;
        end
    end

    %% Define all sizes
    nU=length(u);
    nG=length(G);
    nF=length(F);
    fprintf('    # primal vars = %d, # equal constr = %d, # inequal constr = %d...\n',nU,nG,nF);

    out.u=u;
    out.F=F;
    out.G=G;
    out.nu=nu;
    out.lambda=lambda;

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

    fprintf('    getJ()...');
    t2=clock();
    declareGet(code,f,'getJ__');

    if debugConvergence
        declareGet(code,u,'getU__');
        if nG>0
            declareGet(code,{G,nu},'getGNu__');
        end
        if nF>0
            declareGet(code,{F,lambda},'getFLambda__');
        end
    end

    fprintf('(%.2f sec)\n    1st derivatives...',etime(clock(),t2));
    t2=clock();
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

    %% Compute derivatives
    if trustRegion
        if ~packOptimizationVariables
            error('gradientsWRTvariables not implemented for trustRegion option');
        end
        f_u=gradient(f+.5*addEye2HessianU*norm2(u),u);
        Lf=f+.5*addEye2HessianU*norm2(u);
    else
        if packOptimizationVariables
            f_u=gradient(f,u);
            Lf=f;
        else
            f_u=gradientVector(f_uList,uList);
            Lf=f_uList;
        end
    end
    out.Lf_u=f_u;
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
    end
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

    % gets for derivatives
    declareGet(code,out.Lf_u,'getGrad__');
    declareGet(code,norminf(out.Lf_u),'getNorminf_Grad__');
    if debugConvergence
        declareGet(code,Lf,'getLf__');
        declareGet(code,out.Lf_uu,'getLfuu__');
    end

    %% Barrier variable
    if nF>0
        out.mu=Tvariable('mu__',[],nowarningsamesize,nowarningever);
        muOnes=reshape(out.mu,1);
        muOnes=repmat(muOnes,nF);

        declareSet(code,out.mu,'setMu__');

        % Automatic initialization of lambda
        declareCopy(code,lambda,muOnes./F,'initDualIneq__');

        declareGet(code,{gap,min(F,[],1),min(lambda,[],1)},'getGapMinFMinLambda__');
    else
        gap=Tzeros([]);
        out.mu=Tzeros([]);
        muOnes=Tzeros(0);
    end

    if nG>0
        % Automatic initialization of nu
        % simple
        declareCopy(code,nu,Tones(nG),'initDualEq__');
        % complex
        WW0=[Teye(nU,nU),G_u';G_u,-addEye2HessianEq*Teye(nG,nG)];
        factor_ww0=factor(WW0,[cmexfunction,'_WW0.subscripts'],[cmexfunction,'_WW0.values']);
        if atomicFactorization
            factor_ww0=declareAlias(code,factor_ww0,'factor_ww0',true,nowarningsamesize,nowarningever);
        end
        b0=[F_u'*lambda-f_u;Tzeros(nG)];
        wnu0=factor_ww0\b0;
        declareCopy(code,nu,wnu0(nU+1:end),'initDualEqX__');

        declareGet(code,norminf(G),'getNorminf_G__');
    end

    alphaPrimal=Tvariable('alphaPrimal__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaPrimal,'setAlphaPrimal__');
    alphaDualEq=Tvariable('alphaDualEq__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaDualEq,'setAlphaDualEq__');
    alphaDualIneq=Tvariable('alphaDualIneq__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaDualIneq,'setAlphaDualIneq__');

    fprintf('(%.2f sec)\n    WW...',etime(clock(),t2));
    t2=clock();

    if smallerNewtonMatrix
        %%%%%%%%%%%%%%%%%%
        %% Small matrix %%
        %%%%%%%%%%%%%%%%%%

        LPG=tprod(lambda./F,1,F_u,[1,2]);
        if trustRegion
            WW11=out.Lf_uu+tprod(F_u,[-1,1],LPG,[-1,2],'associate');
            WW=[WW11,G_u';
                G_u,-addEye2HessianEq*Teye([nG,nG])];
        else
            WW11=out.Lf_uu+tprod(F_u,[-1,1],LPG,[-1,2],'associate')+addEye2HessianU*Teye(size(out.Lf_uu));
            WW=[WW11,G_u';
                G_u,-addEye2HessianEq*Teye([nG,nG])];
        end
        out.Hess=WW;
        muF=muOnes./F;         % muF=(mu*Tones(size(F)))./F;

        factor_ww=factor(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if isequal(factor,@ldl)
            out.dHess=ldl_d(factor_ww);
            out.lHess=ldl_l(factor_ww);
            tol=0*1e-8; % minimum eigenvalue to be considered positive -- as we let addEye2Hessian get smaller, this generally also needs to get smaller
            declareGet(code,{sum(heaviside(out.dHess-tol)),...
                             sum(heaviside(-out.dHess-tol))},'getHessInertia__');
            %declareGet(code,{sqrt(norm2(out.Hess-out.lHess*diag(out.dHess)*out.lHess'))},'getFactorError__');
        else
            out.dHess=Tzeros(size(factor_ww,1));
            out.lHess=Tzeros(size(factor_ww));
        end
        if atomicFactorization
            factor_ww=declareAlias(code,factor_ww,'factor_ww',true,nowarningsamesize,nowarningever);
        end

        if skipAffine
            b_s=[-f_u-tprod(G_u,[-1,1],nu,-1)+tprod(F_u,[-1,1],muF,-1);
                 -G];
        else
            %% affine direction
            b_a=[-f_u-tprod(G_u,[-1,1],nu,-1);
                 -G];

            dx_a=factor_ww\b_a;
            dx_a=declareAlias(code,dx_a,'dUNu_a__',false,nowarningsamesize,nowarningever);

            dU_a=dx_a(1:nU);
            newU_a=u+alphaPrimal*dU_a;

            if nF>0
                dLambda_a=-LPG*dU_a-lambda;
                newLambda_a=lambda+alphaDualIneq*dLambda_a;

                newF_a=substitute(F,u,newU_a);

                maxAlphaPrimal_a=clp(F,F_u*dU_a);
                maxAlphaDualIneq_a=clp(lambda,dLambda_a);

                rho=tprod(newF_a,[-1],newLambda_a,[-1])./gap; % rho=(newF_a*newLambda_a)./gap;

                declareGet(code,{maxAlphaPrimal_a,maxAlphaDualIneq_a},'getMaxAlphas_a__');
                declareGet(code,min(newF_a,[],1),'getMinF_a__');
                declareGet(code,rho,'getRho__');
            else
                dLambda_a=Tzeros(nF);
            end
            % Mehrotra correction for search direction
            Mehrotra=(F_u*dU_a).*dLambda_a./F;
            b_s=[-f_u-tprod(G_u,[-1,1],nu,-1)+tprod(F_u,[-1,1],muF,-1)-tprod(F_u,[-1,1],Mehrotra,-1);
                 -G];
        end

        fprintf('(%.2f sec)\n    search directions...',etime(clock(),t2));
        t2=clock();
        %% search direction
        dx_s=factor_ww\b_s;
        dx_s=declareAlias(code,dx_s,'dx_s__',false,nowarningsamesize,nowarningever);

        declareGet(code,{norminf((WW*dx_s-b_s))},'getDirectionError__');

        dU_s=dx_s(1:nU);
        newU_s=u+alphaPrimal*dU_s;
        dNu_s=dx_s(nU+1:end);
        newNu_s=nu+alphaDualEq*dNu_s;
        curvature=dU_s*(WW11*dU_s);
        declareGet(code,curvature,'getCurvature__');

        if nF>0
            if skipAffine
                dLambda_s=muF-LPG*dU_s-lambda;
            else
                dLambda_s=muF-LPG*dU_s-lambda-Mehrotra;
            end
            newLambda_s=lambda+alphaDualIneq*dLambda_s;

            maxAlphaPrimal_s=clp(F,F_u*dU_s);
            maxAlphaDualIneq_s=clp(lambda,dLambda_s);
            declareGet(code,{maxAlphaPrimal_s,maxAlphaDualIneq_s},'getMaxAlphas_s__');

            newF_s=substitute(F,u,newU_s);
            declareGet(code,min(newF_s,[],1),'getMinF_s__');
            if debugConvergence
                declareGet(code,{newF_s,newLambda_s},'getFLambda_s__');
            end
        else
            newLambda_s=Tzeros(nF);
        end

        declareCopy(code,{u,nu,lambda},{newU_s,newNu_s,newLambda_s},'updatePrimalDual__');
    else
        %%%%%%%%%%%%%%%%%%
        %% Large matrix %%
        %%%%%%%%%%%%%%%%%%

        if trustRegion
            WW11=out.Lf_uu;
            WW=[WW11,G_u',-F_u';
                G_u,-addEye2HessianEq*Teye([nG,nG]),Tzeros([nG,nF]);
                -F_u,Tzeros([nF,nG]),-diag(F./lambda)];
        else
            WW11=out.Lf_uu+addEye2HessianU*Teye(size(out.Lf_uu));
            WW=[WW11,G_u',-F_u';
                G_u,-addEye2HessianEq*Teye([nG,nG]),Tzeros([nG,nF]);
                -F_u,Tzeros([nF,nG]),-diag(F./lambda)];
        end
        out.Hess=WW;

        factor_ww=factor(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if isequal(factor,@ldl)
            out.dHess=ldl_d(factor_ww);
            out.lHess=ldl_l(factor_ww);
            tol=0*1e-8; % minimum eigenvalue to be considered positive -- as we let addEye2Hessian get smaller, this generally also needs to get smaller
            declareGet(code,{sum(heaviside(out.dHess-tol)),...
                             sum(heaviside(-out.dHess-tol))},'getHessInertia__');
            %declareGet(code,{sqrt(norm2(out.Hess-out.lHess*diag(out.dHess)*out.lHess'))},'getFactorError__');
        else
            out.dHess=Tzeros(size(factor_ww,1));
            out.lHess=Tzeros(size(factor_ww));
        end
        if atomicFactorization
            factor_ww=declareAlias(code,factor_ww,'factor_ww',true,nowarningsamesize,nowarningever);
        end

        if skipAffine
            b_s=[-out.Lf_u;
                 -G;
                 F-muOnes./lambda];
        else
            %% affine direction
            b_a=[-out.Lf_u;
                 -G;
                 F];

            dx_a=factor_ww\b_a;
            dx_a=declareAlias(code,dx_a,'dUNu_a__',false,nowarningsamesize,nowarningever);

            dU_a=dx_a(1:nU);
            newU_a=u+alphaPrimal*dU_a;

            if nF>0
                dLambda_a=dx_a(nU+nG+1:nU+nG+nF);
                newLambda_a=lambda+alphaDualIneq*dLambda_a;

                newF_a=substitute(F,u,newU_a);
                maxAlphaPrimal_a=clp(F,F_u*dU_a);
                maxAlphaDualIneq_a=clp(lambda,dLambda_a);

                rho=tprod(newF_a,[-1],newLambda_a,[-1])./gap; % rho=(newF_a*newLambda_a)./gap;

                declareGet(code,{maxAlphaPrimal_a,maxAlphaDualIneq_a},'getMaxAlphas_a__');

                declareGet(code,min(newF_a,[],1),'getMinF_a__');
                declareGet(code,rho,'getRho__');
            else
                dLambda_a=Tzeros(nF);
            end

            % Mehrotra correction for search direction
            b_s=[-out.Lf_u;
                 -G;
                 F+(F_u*dU_a).*dLambda_a./lambda-muOnes./lambda];
        end

        %% search direction
        dx_s=factor_ww\b_s;
        dx_s=declareAlias(code,dx_s,'dx_s__',false,nowarningsamesize,nowarningever);

        declareGet(code,{norminf((WW*dx_s-b_s))},'getDirectionError__');

        dU_s=dx_s(1:nU);
        newU_s=u+alphaPrimal*dU_s;
        dNu_s=dx_s(nU+1:nU+nG);
        newNu_s=nu+alphaDualEq*dNu_s;

        curvature=dU_s*(WW11*dU_s);
        declareGet(code,curvature,'getCurvature__');

        if nF>0
            dLambda_s=dx_s(nU+nG+1:nU+nG+nF);
            newLambda_s=lambda+alphaDualIneq*dLambda_s;

            maxAlphaPrimal_s=clp(F,F_u*dU_s);
            maxAlphaDualIneq_s=clp(lambda,dLambda_s);
            declareGet(code,{maxAlphaPrimal_s,maxAlphaDualIneq_s},'getMaxAlphas_s__');
            newF_s=substitute(F,u,newU_s);
            declareGet(code,min(newF_s,[],1),'getMinF_s__');
            if debugConvergence
                declareGet(code,{newF_s,newLambda_s},'getFLambda_s__');
            end
        else
            newLambda_s=Tzeros(nF);
        end
        declareCopy(code,{u,nu,lambda},{newU_s,newNu_s,newLambda_s},'updatePrimalDual__');
    end % smallerNewtonMatrix

    out.b_s=b_s;
    out.dx_s=dx_s;


    if debugConvergence;
        declareGet(code,{dU_s,dNu_s,dLambda_s},'getD__');
        declareGet(code,out.Hess,'getHess__');
        declareGet(code,full(out.dHess),'getdHess__');
    end

    if any(isSensitivity)
        if smallerNewtonMatrix
            error('ipmPD_CS: computation of sensitivity not implemented for smallerNewtonMatrix\n');
        end
        out.DfDu1=out.Lf_u(isSensitivity);
        out.D2fDu1=WW([isSensitivity;false(nF+nG,1)],[isSensitivity;false(nF+nG,1)]);
        if any(~isSensitivity)
            out.Hess1 =  WW([~isSensitivity;true(nG+nF,1)],[~isSensitivity;true(nG+nF,1)]);
            out.B1    = -WW([~isSensitivity;true(nG+nF,1)],[isSensitivity;false(nG+nF,1)]);
            if atomicFactorization
                %factor_Hess1=declareAlias(code,factor_Hess1,'factor_Hess1',true,nowarningsamesize,nowarningever);
                % atomic mldivide(LU,b) only implemented for 1-d vector
                out.Du1=Tzeros(size(out.B1));
            else
                factor_Hess1=factor(out.Hess1);
                out.Du1=factor_Hess1\out.B1;
            end
            out.D2fDu1=out.D2fDu1+WW([isSensitivity;false(nF+nG,1)],[~isSensitivity;true(nG+nF,1)])*out.Du1;
        else
            out.Du1=Tconstant([]);
        end
    else
        out.Du1=Tconstant([]);
        out.DfDu1=Tconstant([]);
        out.D2fDu1=Tconstant([]);
    end

    % declareGet(code,full(WW),'getWW__');
    % declareGet(code,u,'getU__');
    % declareGet(code,lambda,'getLambda__');
    % declareGet(code,b_s,'getb_s__');
    % declareGet(code,dx_s,'getDx_s__');

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
    fprintf('  done ipmPD_CS symbolic computations (%.3f sec)\n',etime(clock(),t1));

    %    profile off;profile viewer

end
