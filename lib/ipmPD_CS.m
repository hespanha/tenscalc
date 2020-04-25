function [Hess_,dHess_,Lf_u,mu,Du1__,DfDu1__,D2fDu1__]=ipmPD_CS(code,f,u,lambda,nu,F,G,isSensitivity,...
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

%profile on

    nowarningsamesize=true;
    nowarningever=true;

    szHess_=TcheckVariable('Hess_');

    if smallerNewtonMatrix
        fprintf('\n  Starting ipmPD_CS symbolic computations (smallNewtonMatrix)...\n');
    else
        fprintf('\n  Starting ipmPD_CS symbolic computations (largeNewtonMatrix)...\n');
    end
    t1=clock();
    
    %% Define all sizes
    nU=length(u);
    nG=length(G);
    nF=length(F);
    fprintf('    # primal vars = %d, # equal constr = %d, # inequal constr = %d...\n',nU,nG,nF);
    
    
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
    Lf_u=f_u;
    
    if nF>0
        mu=Tvariable('mu__',[],nowarningsamesize,nowarningever);
        %muOnes=mu*Tones(nF);
        muOnes=reshape(mu,1);
        muOnes=muOnes(ones(nF,1));
        
        declareSet(code,mu,'setMu__');
        
        F_u=gradient(F,u);
        gap=tprod(lambda,-1,F,-1);                    % gap=lambda*F;
        Lf=Lf-gap;                                    % Lf=Lf-gap;
        Lf_u=Lf_u-tprod(F_u,[-1,1],lambda,-1);        % Lf_u=Lf_u-F_u'*lambda;
        
        % Automatic initialization of lambda
        declareCopy(code,lambda,muOnes./F,'initDualIneq__');
        
        declareGet(code,{gap,min(F,[],1),min(lambda,[],1)},'getGapMinFMinLambda__');
    else
        F=Tzeros(0);
        F_u=Tzeros([0,nU]);
        gap=Tzeros([]);
        mu=Tzeros([]);
        muOnes=Tzeros(0);
    end

    if nG>0
        G_u=gradient(G,u);
        Lf=Lf+tprod(nu,-1,G,-1);                      % Lf=Lf+nu*G;
        Lf_u=Lf_u+tprod(G_u,[-1,1],nu,-1);            % Lf_u=Lf_u+G_u'*nu;

        % Automatic initialization of nu
        declareCopy(code,nu,Tones(nG),'initDualEq__');

        declareGet(code,norminf(G),'getNorminf_G__');
    else
        G=Tzeros(0);
        G_u=Tzeros([0,nU]);
    end
    
    declareGet(code,Lf_u,'getGrad__');
    declareGet(code,norminf(Lf_u),'getNorminf_Grad__');

    fprintf('(%.2f sec)\n    2nd derivatives...',etime(clock(),t2));
    t2=clock();

    Lf_uu=gradient(Lf_u,u);
    
    if debugConvergence
        declareGet(code,Lf,'getLf__');
        declareGet(code,full(Lf_uu),'getLfuu__');
    end
    
    alphaPrimal=Tvariable('alphaPrimal__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaPrimal,'setAlphaPrimal__');
    alphaDualEq=Tvariable('alphaDualEq__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaDualEq,'setAlphaDualEq__');
    alphaDualIneq=Tvariable('alphaDualIneq__',[],nowarningsamesize,nowarningever);
    declareSet(code,alphaDualIneq,'setAlphaDualIneq__');

    if useLDL
        if atomicFactorization
            factor=@lu_sym;
        else
            factor=@ldl;
        end        
    else
        factor=@lu;
    end
    
    fprintf('(%.2f sec)\n    WW...',etime(clock(),t2));
    t2=clock();
    if smallerNewtonMatrix
        %%%%%%%%%%%%%%%%%%
        %% Small matrix %%
        %%%%%%%%%%%%%%%%%%

        LPG=tprod(lambda./F,1,F_u,[1,2]);
        Hess_=[Lf_uu+tprod(F_u,[-1,1],LPG,[-1,2],'associate'),G_u';
               G_u,Tzeros([nG,nG])];
        WW=  [Lf_uu+tprod(F_u,[-1,1],LPG,[-1,2],'associate')+tprod(addEye2Hessian,[],Teye(size(Lf_uu)),[1,2]),G_u';
              G_u,tprod(-addEye2Hessian,[],Teye([nG,nG]),[1,2])];
        muF=muOnes./F;         % muF=(mu*Tones(size(F)))./F;
        
        factor_ww=factor(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if isequal(factor,@ldl)
            dHess_=ldl_d(factor_ww);
        else
            dHess_=Tzeros(size(factor_ww,1));
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

        Hess_=[Lf_uu,G_u',-F_u';
               G_u,Tzeros([nG,nG+nF]);
               -F_u,Tzeros([nF,nG]),-diag(F./lambda)];
        WW=[Lf_uu+tprod(addEye2Hessian,[],Teye(size(Lf_uu)),[1,2]),G_u',-F_u';
            G_u,-tprod(addEye2Hessian,[],Teye([nG,nG]),[1,2]),Tzeros([nG,nF]);
            -F_u,Tzeros([nF,nG]),-diag(F./lambda)-tprod(addEye2Hessian,[],Teye([nF,nF]),[1,2])];

        factor_ww=factor(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if isequal(factor,@ldl)
            dHess_=ldl_d(factor_ww);
        else
            dHess_=Tzeros(size(factor_ww,1));
        end
        if atomicFactorization
            factor_ww=declareAlias(code,factor_ww,'factor_ww',true,nowarningsamesize,nowarningever);
        end

        if skipAffine
            b_s=[-Lf_u;
                 -G;
                 F-muOnes./lambda];
        else
            %% affine direction
            b_a=[-Lf_u;
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
            b_s=[-Lf_u;
                 -G;
                 F+(F_u*dU_a).*dLambda_a./lambda-muOnes./lambda];
        end
        
        %% search direction
        dx_s=factor_ww\b_s;
        dx_s=declareAlias(code,dx_s,'dx_s__',false,nowarningsamesize,nowarningever);
        
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
    
    if debugConvergence;
        declareGet(code,Hess_,'getHess__');
    end

    if nargout>2 && any(isSensitivity)
        if smallerNewtonMatrix
            error('ipmPD_CS: computation of sensitivity not implemented for smallerNewtonMatrix\n');
        end
        DfDu1__=Lf_u(isSensitivity);
        D2fDu1__=WW([isSensitivity;false(nF+nG,1)],[isSensitivity;false(nF+nG,1)]);        
        if any(~isSensitivity)
            Hess1=WW([~isSensitivity;true(nG+nF,1)],[~isSensitivity;true(nG+nF,1)]);
            B1=-WW([~isSensitivity;true(nG+nF,1)],[isSensitivity;false(nG+nF,1)]);
            if atomicFactorization
                %factor_Hess1=declareAlias(code,factor_Hess1,'factor_Hess1',true,nowarningsamesize,nowarningever);
                % atomic mldivide(LU,b) only implemented for 1-d vector
                Du1__=Tzeros(size(B1));
            else
                factor_Hess1=factor(Hess1);
                Du1__=factor_Hess1\B1;
            end
            D2fDu1__=D2fDu1__+WW([isSensitivity;false(nF+nG,1)],[~isSensitivity;true(nG+nF,1)])*Du1__;        
        else
            Du1__=Tconstant([]);
        end
    else
        Du1__=Tconstant([]);
        DfDu1__=Tconstant([]);
        D2fDu1__=Tconstant([]);
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

    if ~isempty(szHess_) && ~myisequal(szHess_,size(Hess_))
        warning('\nvariable: ''Hess_'' already exists with the wrong size [%d,%d], should be [%d,%d]\n',...
                szHess_(1),szHess_(2),size(Hess_,1),size(Hess_,2));
    end

    fprintf('(%.2f sec)\n    ',etime(clock(),t2));
    fprintf('  done ipmPD_CS symbolic computations (%.3f sec)\n',etime(clock(),t1));
    
    %profile off
    %profile viewer
    
end

