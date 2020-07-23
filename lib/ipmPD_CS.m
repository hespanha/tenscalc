function out=ipmPD_CS(code,f,u,lambda,nu,F,G,isSensitivity,...
                      smallerNewtonMatrix,addEye2Hessian,skipAffine,...
                      scaleInequalities,scaleCost,scaleEqualities,...
                      useLDL,atomicFactorization,...
                      cmexfunction,allowSave,debugConvergence)
% See ../doc/ipm.tex for an explanation of the formulas used here
%
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

%profile clear;profile on

    nowarningsamesize=true;
    nowarningever=true;

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
        else
            factor=@ldl;
        end        
    else
        factor=@lu;
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
    end
    if scaleCost>0
        scale4Cost=Tvariable('scale4Cost__',[]);
        declareCopy(code,scale4Cost,abs(scaleCost/f),'scaleCost__');
        f=scale4Cost*f;
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
    
    fprintf('(%.2f sec)\n    1st derivates...',etime(clock(),t2));
    t2=clock();
    f_u=gradient(f,u);
    Lf=f;
    out.Lf_u=f_u;
    
    if addEye2Hessian
        addEye2Hessian1=Tvariable('addEye2Hessian1__',[],nowarningsamesize,nowarningever);
        addEye2Hessian2=Tvariable('addEye2Hessian2__',[],nowarningsamesize,nowarningever);
        declareSet(code,addEye2Hessian1,'setAddEye2Hessian1__');
        declareSet(code,addEye2Hessian2,'setAddEye2Hessian2__');
    else
        addEye2Hessian1=Tzeros([]);
        addEye2Hessian2=Tzeros([]);
    end
    
    if nF>0
        out.mu=Tvariable('mu__',[],nowarningsamesize,nowarningever);
        %muOnes=out.mu*Tones(nF);
        muOnes=reshape(out.mu,1);
        muOnes=muOnes(ones(nF,1));
        
        declareSet(code,out.mu,'setMu__');
        
        F_u=gradient(F,u);
        gap=tprod(lambda,-1,F,-1);                    % gap=lambda*F;
        Lf=Lf-gap;                                    % Lf=Lf-gap;
        out.Lf_u=out.Lf_u-tprod(F_u,[-1,1],lambda,-1);        % Lf_u=Lf_u-F_u'*lambda;
        
        % Automatic initialization of lambda
        declareCopy(code,lambda,muOnes./F,'initDualIneq__');
        
        declareGet(code,{gap,min(F,[],1),min(lambda,[],1)},'getGapMinFMinLambda__');
    else
        F=Tzeros(0);
        F_u=Tzeros([0,nU]);
        gap=Tzeros([]);
        out.mu=Tzeros([]);
        muOnes=Tzeros(0);
    end

    if nG>0
        G_u=gradient(G,u);
        Lf=Lf+tprod(nu,-1,G,-1);                      % Lf=Lf+nu*G;
        out.Lf_u=out.Lf_u+tprod(G_u,[-1,1],nu,-1);            % Lf_u=Lf_u+G_u'*nu;

        % Automatic initialization of nu
        % simple
        declareCopy(code,nu,Tones(nG),'initDualEq__');
        % complex
        WW0=[Teye(nU,nU),G_u';G_u,-addEye2Hessian2*Teye(nG,nG)];
        factor_ww0=factor(WW0,[cmexfunction,'_WW0.subscripts'],[cmexfunction,'_WW0.values']);
        if atomicFactorization
            factor_ww0=declareAlias(code,factor_ww0,'factor_ww0',true,nowarningsamesize,nowarningever);
        end
        b0=[F_u'*lambda-f_u;Tzeros(nG)];
        wnu0=factor_ww0\b0;
        declareCopy(code,nu,wnu0(nU+1:end),'initDualEqX__');
        
        declareGet(code,norminf(G),'getNorminf_G__');
    else
        G=Tzeros(0);
        G_u=Tzeros([0,nU]);
    end
    
    declareGet(code,out.Lf_u,'getGrad__');
    declareGet(code,norminf(out.Lf_u),'getNorminf_Grad__');

    fprintf('(%.2f sec)\n    2nd derivatives...',etime(clock(),t2));
    t2=clock();

    out.Lf_uu=gradient(out.Lf_u,u);
    
    if debugConvergence
        declareGet(code,Lf,'getLf__');
        declareGet(code,out.Lf_uu,'getLfuu__');
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
        WW=  [out.Lf_uu+tprod(F_u,[-1,1],LPG,[-1,2],'associate')+addEye2Hessian1*Teye(size(out.Lf_uu)),G_u';
              G_u,-addEye2Hessian2*Teye([nG,nG])];
        out.Hess=WW;
        muF=muOnes./F;         % muF=(mu*Tones(size(F)))./F;
        
        factor_ww=factor(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if isequal(factor,@ldl)
            out.dHess=ldl_d(factor_ww);
            out.lHess=ldl_l(factor_ww);
            tol=0*1e-8; % minimum eigenvalue to be considered positive -- as we let addEye2Hessian get smaller, this generally also needs to get smaller
            declareGet(code,{sum(heaviside(out.dHess-tol)),sum(heaviside(-out.dHess-tol))},'getHessInertia__');
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

        WW=[out.Lf_uu+addEye2Hessian1*Teye(size(out.Lf_uu)),G_u',-F_u';
            G_u,-addEye2Hessian2*Teye([nG,nG]),Tzeros([nG,nF]);
            -F_u,Tzeros([nF,nG]),-diag(F./lambda)];
        out.Hess=WW;
        
        factor_ww=factor(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if isequal(factor,@ldl)
            out.dHess=ldl_d(factor_ww);
            out.lHess=ldl_l(factor_ww);
            tol=0*1e-8; % minimum eigenvalue to be considered positive -- as we let addEye2Hessian get smaller, this generally also needs to get smaller
            declareGet(code,{sum(heaviside(out.dHess-tol)),sum(heaviside(-out.dHess-tol))},'getHessInertia__');
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

