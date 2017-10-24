function ipmPDeqlat_CS(code,f,g,u,d,x,P1lambda,P1nu,P1xnu,P2lambda,P2nu,P2xnu,...
                       Fu,Gu,Fd,Gd,H,...
                       smallerNewtonMatrix,addEye2Hessian,skipAffine,...
                       atomicFactorization,...
                       cmexfunction,allowSave,debugConvergence,profiling)
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% See ../doc/ipm.tex for an explanation of the formulas used here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%profile on
    
    if smallerNewtonMatrix
        fprintf('\n  Starting ipmPDeqlat_CS symbolic computations (smallNewtonMatrix)...\n');
    else
        fprintf('\n  Starting ipmPDeqlat_CS symbolic computations (largeNewtonMatrix)...\n');
    end

    t1=clock();
    
    %% Define all sizes
    nU=length(u);
    nD=length(d);
    nX=length(x);
    nZ=nU+nD+nX;
    nGu=length(Gu);
    nGd=length(Gd);
    nH=length(H);
    nG=nGu+nGd+nH;
    nFu=length(Fu);
    nFd=length(Fd);
    nF=nFu+nFd;
    
    fprintf('    getfg()...');
    t2=clock();
    declareGet(code,{f,g},'getfg__');

    %% Stack variables
    F=[Fu;Fd];
    lambda=[P1lambda;P2lambda];
    G=[];
    nu=[];
    if nGu>0
        G=[G;Gu];
        nu=[nu;P1nu];
    end
    if nH>0
        nu=[nu;P1xnu];
    end
    if nGd>0
        G=[G;Gd];
        nu=[nu;P2nu];
    end
    if nH>0
        G=[G;H];
        nu=[nu;P2xnu];
    end
    gap=lambda*F;
    
    if debugConvergence
        declareGet(code,{u,d},'getUD__');
        if nG>0
            declareGet(code,{G,nu},'getGNu__');
        end
        if nF>0
            declareGet(code,{F,lambda},'getFLambda__');
        end
    end 
    
    %% Declare gets for exit condition and output
    if nF>0
        mu=Tvariable('mu__',[]);
        %muOnes=mu*Tones(nF);
        muOnes=reshape(mu,1);
        muOnes=muOnes(ones(nF,1));

        declareSet(code,mu,'setMu__');
        declareGet(code,{gap,min(F,1),min(lambda,1)},'getGapMinFMinLambda__');
    end
    if nG>0
        declareGet(code,norminf(G),'getNorminf_G__');
    end
    
    %% Automatic initialization of the lambda variable
    if 0
        declareCopy(code,{nu,lambda},{Tones(nG+nH),mu*Tones(nF)./F},'initDual__');
    else
        dst={};
        src={};
        if nGu>0
            dst{end+1}=P1nu;
            src{end+1}=Tones(nGu);
        end
        if nGd>0
            dst{end+1}=P2nu;
            src{end+1}=Tones(nGd);
        end
        if nH>0
            dst{end+1}=P1xnu;
            src{end+1}=Tones(nH);
            dst{end+1}=P2xnu;
            src{end+1}=Tones(nH);
        end
        if ~isempty(dst)
            declareCopy(code,dst,src,'initDualEq__');
        end
        if nF>0
            dst={};
            src={};
            if nFu>0
                dst{end+1}=P1lambda;
                src{end+1}=mu*Tones(nFu)./Fu;
            end
            if nFd>0
                dst{end+1}=P2lambda;
                src{end+1}=mu*Tones(nFd)./Fd;
            end
            declareCopy(code,dst,src,'initDualIneq__');
        end
    end
    
    %% Construct the "Hessian" matrix
    Lf=f;
    Lg=g;
    if nFu>0
        Lf=Lf-P1lambda*Fu;
    end
    if nFd>0
        Lg=Lg-P2lambda*Fd;
    end
    if nGu>0
        Lf=Lf+P1nu*Gu;
    end
    if nGd>0
        Lg=Lg+P2nu*Gd;
    end
    if nH>0
        Lf=Lf+P1xnu*H;
        Lg=Lg+P2xnu*H;
    end
    
    alphaPrimal=Tvariable('alphaPrimal__',[]);
    declareSet(code,alphaPrimal,'setAlphaPrimal__');
    alphaDualEq=Tvariable('alphaDualEq__',[]);
    declareSet(code,alphaDualEq,'setAlphaDualEq__');
    alphaDualIneq=Tvariable('alphaDualIneq__',[]);
    declareSet(code,alphaDualIneq,'setAlphaDualIneq__');

    fprintf('(%.2f sec)\n    1st derivates...',etime(clock(),t2));
    t2=clock();
    Lf_u=gradient(Lf,u);
    Lg_d=gradient(Lg,d);
    if nX>0
        Lf_x=gradient(Lf,x);
        Lg_x=gradient(Lg,x);
        % for exit condition
        declareGet(code,norminf(Lf_u)+norminf(Lf_x)+norminf(Lg_d)+norminf(Lg_x),'getNorminf_Grad__'); 
        % declareGet(code,norminf([gradient(f+P1nu*Gu+P1xnu*H,u);
        %                       gradient(f+P1nu*Gu+P1xnu*H,x);
        %                       gradient(g+P2nu*Gd+P2xnu*H,d);
        %                       gradient(g+P2nu*Gd+P2xnu*H,x);]),'getNorminf_Grad__'); 
    else
        Lf_x=Tzeros(0);
        Lg_x=Tzeros(0);
        % for exit condition
        declareGet(code,norminf(Lf_u)+norminf(Lg_d),'getNorminf_Grad__'); 
    end

    fprintf('(%.2f sec)\n    2nd derivatives...',etime(clock(),t2));
    t2=clock();
    % Derivatives needed to compute WW
    if nX>0
        Lf_uz=gradientVector(Lf_u,{u,d,x});
        Lf_xz=gradientVector(Lf_x,{u,d,x});
        Lg_dz=gradientVector(Lg_d,{u,d,x});
        Lg_xz=gradientVector(Lg_x,{u,d,x});
        
        Lf_ul=gradientVector(Lf_u,{P1lambda,P2lambda});
        Lf_xl=gradientVector(Lf_x,{P1lambda,P2lambda});
        Lg_dl=gradientVector(Lg_d,{P1lambda,P2lambda});
        Lg_xl=gradientVector(Lg_x,{P1lambda,P2lambda});
        
        Lf_un=gradientVector(Lf_u,{P1nu,P1xnu,P2nu,P2xnu});
        Lf_xn=gradientVector(Lf_x,{P1nu,P1xnu,P2nu,P2xnu});
        Lg_dn=gradientVector(Lg_d,{P1nu,P1xnu,P2nu,P2xnu});
        Lg_xn=gradientVector(Lg_x,{P1nu,P1xnu,P2nu,P2xnu});
        
        F_z=gradientVector(F,{u,d,x});
        G_z=gradientVector(G,{u,d,x});
    else
        Lf_uz=gradientVector(Lf_u,{u,d});
        Lg_dz=gradientVector(Lg_d,{u,d});
        
        Lf_ul=gradientVector(Lf_u,{P1lambda,P2lambda});
        Lg_dl=gradientVector(Lg_d,{P1lambda,P2lambda});
        
        Lf_un=gradientVector(Lf_u,{P1nu,P1xnu,P2nu,P2xnu});
        Lg_dn=gradientVector(Lg_d,{P1nu,P1xnu,P2nu,P2xnu});
        
        F_z=gradientVector(F,{u,d});
        G_z=gradientVector(G,{u,d});
    end

    if nF>0
        LFF=tprod(lambda./F,1,F_z,[1,2]); %LFF=diag(lambda./F)*F_z;

        fprintf('(%.2f sec)\n    adding LFF...',etime(clock(),t2));
        t2=clock();
        LFF=declareAlias(code,LFF,'LFF__');
    end
    
    fprintf('(%.2f sec)\n    WW & b...',etime(clock(),t2));
    t2=clock();
    if smallerNewtonMatrix
        %%%%%%%%%%%%%%%%%%
        %% Small matrix %%
        %%%%%%%%%%%%%%%%%%

        if nX>0
            if nF>0
                WW=[Lf_uz-Lf_ul*LFF,Lf_un;
                    Lg_dz-Lg_dl*LFF,Lg_dn;
                    Lf_xz-Lf_xl*LFF,Lf_xn;
                    Lg_xz-Lg_xl*LFF,Lg_xn;
                    G_z,Tzeros([nG,nG+nH])];
            else
                WW=[Lf_uz,Lf_un;
                    Lg_dz,Lg_dn;
                    Lf_xz,Lf_xn;
                    Lg_xz,Lg_xn;
                    G_z,Tzeros([nG,nG+nH])];
            end
            %WW=WW+addEye2Hessian*Teye([nZ+nH+nG,nZ+nH+nG]);
            
            ff=f;
            gg=g;
            if nGu>0
                ff=ff+P1nu*Gu;
            end
            if nGd>0
                gg=gg+P2nu*Gd;
            end
            if nH>0
                ff=ff+P1xnu*H;
                gg=gg+P2xnu*H;
            end
            
            b_a=-[gradient(ff,u);
                  gradient(gg,d);
                  gradient(ff,x);
                  gradient(gg,x);
                  G];
            b_s=b_a;
            if nF>0
                b_s=b_s-[Lf_ul*(muOnes./F);
                         Lg_dl*(muOnes./F);
                         Lf_xl*(muOnes./F);
                         Lg_xl*(muOnes./F);
                         Tzeros(nG)];
            end
        else
            if nF>0
                WW=[[Lf_uz-Lf_ul*LFF;
                     Lg_dz-Lg_dl*LFF]+addEye2Hessian*Teye((nZ)*[1,1]),...
                    [Lf_un;
                     Lg_dn];
                    G_z,-addEye2Hessian*Teye((nG)*[1,1])];
                %WW=[Lf_uz-Lf_ul*LFF,Lf_un;
                %    Lg_dz-Lg_dl*LFF,Lg_dn;
                %    G_z,Tzeros([nG,nG+nH])];
                %WW=WW+addEye2Hessian*Teye([nZ+nH+nG,nZ+nH+nG]);
            else
                WW=[[Lf_uz;
                     Lg_dz]+addEye2Hessian*Teye((nZ)*[1,1]),...
                    [Lf_un;
                     Lg_dn];
                    G_z,-addEye2Hessian*Teye((nG)*[1,1])];
            end
            
            ff=f;
            gg=g;
            if nGu>0
                ff=ff+P1nu*Gu;
            end
            if nGd>0
                gg=gg+P2nu*Gd;
            end
            % if nH>0
            %     ff=ff+P1xnu*H;
            %     gg=gg+P2xnu*H;
            % end
            
            b_a=-[gradient(ff,u);
                  gradient(gg,d);
                  G];
            
            b_s=b_a;
            if nF>0
                b_s=b_s-[Lf_ul*(muOnes./F);
                         Lg_dl*(muOnes./F);
                         Tzeros(nG)];
            end
        end
        
        fprintf('(%.2f sec)\n    adding WW...',etime(clock(),t2));
        t2=clock();

        WW=declareAlias(code,WW,'WW__');
        
        if debugConvergence
            declareGet(code,{full([u;d]),full(nu),full(lambda)},'getZNL__');
            declareGet(code,{full(F),full(G)},'getFG__');
        end

        lu_ww=lu(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if atomicFactorization
            lu_ww=declareAlias(code,lu_ww,'lu_ww',true);
        end
        
        if ~skipAffine
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% affine scaling direction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('(%.2f sec)\n    affine direction...',etime(clock(),t2));
            t2=clock();
            
            dZNu_a=lu_ww\b_a;
            
            fprintf('(%.2f sec)\n    updated vectores...',etime(clock(),t2));
            t2=clock();
            
            dZNu_a=declareAlias(code,dZNu_a,'dZNu_a__');		    
            dU_a=dZNu_a(1:nU);
            dD_a=dZNu_a(nU+1:nU+nD);
            dZ_a=dZNu_a(1:nZ);
            newU_a=u+alphaPrimal*dU_a;
            newD_a=d+alphaPrimal*dD_a;
            if nX>0
                dX_a=dZNu_a(nU+nD+1:nZ);
                newX_a=x+alphaPrimal*dX_a;
            end
            
            if nF>0
                dLambda_a=-lambda-LFF*dZ_a;
                newLambda_a=lambda+alphaDualIneq*dLambda_a;
                
                newF_a=substitute(F,u,newU_a);
                newF_a=substitute(newF_a,d,newD_a);
                if nX>0
                    newF_a=substitute(newF_a,x,newX_a);
                end
                
                maxAlphaPrimal_a=clp(F,F_z*dZ_a);
                maxAlphaDual_a=clp(lambda,dLambda_a);
                
                % Mehrotra correction
                Mehrotra=(F_z*dZ_a).*dLambda_a./F;
                if nX>0
                    b_s=b_s+[Lf_ul*Mehrotra;
                             Lg_dl*Mehrotra;
                             Lf_xl*Mehrotra;
                             Lg_xl*Mehrotra;
                             Tzeros(nG)];
                else
                    b_s=b_s+[Lf_ul*Mehrotra;
                             Lg_dl*Mehrotra;
                             Tzeros(nG)];
                end
                
                rho=(newF_a*newLambda_a)./gap;
                
                declareGet(code,{maxAlphaPrimal_a,maxAlphaDualIneq_a},'getMaxAlphas_a__');
                declareGet(code,min(newF_a,1),'getMinF_a__');
                declareGet(code,rho,'getRho__');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% combined search direction  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('(%.2f sec)\n    combined direction...',etime(clock(),t2));
        t2=clock();
        
        dZNu_s=lu_ww\b_s;
        
        
        fprintf('(%.2f sec)\n    add direction...',etime(clock(),t2));
        t2=clock();
        
        dZNu_s=declareAlias(code,dZNu_s,'dZNu_s__');		    
        
        fprintf('(%.2f sec)\n    increment primal...',etime(clock(),t2));
        t2=clock();
        
        dU_s=dZNu_s(1:nU);
        dD_s=dZNu_s(nU+1:nU+nD);
        dZ_s=dZNu_s(1:nZ);
        
        fprintf('(%.2f sec)\n    increment dual...',etime(clock(),t2));
        t2=clock();
        
        if nX>0
            dX_s=dZNu_s(nU+nD+1:nZ);
        end
        
        dNu_s=dZNu_s(nZ+1:end);
        dNuU_s=dZNu_s(nZ+1:nZ+nGu);
        dNuD_s=dZNu_s(nZ+nGu+nH+1:nZ+nGu+nH+nGd);
        dNuUx_s=dZNu_s(nZ+nGu+1:nZ+nGu+nH);
        dNuDx_s=dZNu_s(nZ+nGu+nH+nGd+1:end);
        
        if nF>0
            dLambda_s=-lambda-LFF*dZ_s+muOnes./F;
        else
            dLambda_s=[];
        end
        
        if ~skipAffine && nF>0
            % Mehrotra correction
            dLambda_s=dLambda_s-Mehrotra;
        end
        
        fprintf('(%.2f sec)\n    new vectores...',etime(clock(),t2));
        t2=clock();
        
        newU_s=u+alphaPrimal*dU_s;
        newD_s=d+alphaPrimal*dD_s;
        if nX>0
            newX_s=x+alphaPrimal*dX_s;
        end
        
        if nF>0
            maxAlphaPrimal_s=clp(F,F_z*dZ_s);
            maxAlphaDual_s=clp(lambda,dLambda_s);
            newLambda_s=lambda+alphaDualIneq*dLambda_s;
            newF_s=substitute(F,u,newU_s);
            newF_s=substitute(newF_s,d,newD_s);
            if nX>0
                newF_s=substitute(newF_s,x,newX_s);
            end
            declareGet(code,{maxAlphaPrimal_s,maxAlphaDual_s},'getMaxAlphas_s__');
            declareGet(code,min(newF_s,1),'getMinF_s__');
            if debugConvergence
                declareGet(code,{newF_s,newLambda_s},'getFLambda_s__');
            end
        end
        
    else % smallerNewtonMatrix
        %%%%%%%%%%%%%%%%%%
        %% Large matrix %%
        %%%%%%%%%%%%%%%%%%

        if nX>0
            if nF>0
                WW=[Lf_uz,Lf_un,Lf_ul;
                    Lg_dz,Lg_dn,Lg_dl;
                    Lf_xz,Lf_xn,Lf_xl;
                    Lg_xz,Lg_xn,Lg_xl;
                    G_z,Tzeros([nG,nG+nH+nF]);
                    F_z,Tzeros([nF,nG+nH]),diag(F./lambda)];
            else
                WW=[Lf_uz,Lf_un;
                    Lg_dz,Lg_dn;
                    Lf_xz,Lf_xn;
                    Lg_xz,Lg_xn;
                    G_z,Tzeros([nG,nG+nH])];
            end
            
            if nF>0
                b_a=[-Lf_u;
                     -Lg_d;
                     -Lf_x;
                     -Lg_x;
                     -G;
                     -F];
                b_s=[-Lf_u;
                     -Lg_d;
                     -Lf_x;
                     -Lg_x;
                     -G;
                     -F+muOnes./lambda];
            else
                b_a=[-Lf_u;
                     -Lg_d;
                     -Lf_x;
                     -Lg_x;
                     -G];
                b_s=b_a;
            end
        else % nX>0
            if nF>0
                WW=[Lf_uz,Lf_un,Lf_ul;
                    Lg_dz,Lg_dn,Lg_dl;
                    G_z,Tzeros([nG,nG+nH+nF]);
                    F_z,Tzeros([nF,nG+nH]),diag(F./lambda)];
            else
                WW=[Lf_uz,Lf_un;
                    Lg_dz,Lg_dn;
                    G_z,Tzeros([nG,nG+nH])];
            end

            if nF>0
                b_a=[-Lf_u;
                     -Lg_d;
                     -G;
                     -F];
                b_s=[-Lf_u;
                     -Lg_d;
                     -G;
                     -F+muOnes./lambda];
            else
                b_a=[-Lf_u;
                     -Lg_d;
                     -G];
                b_s=b_a;
            end
        end
        
        fprintf('(%.2f sec)\n    adding WW...',etime(clock(),t2));
        t2=clock();

        WW=declareAlias(code,WW,'WW__');
        
        if debugConvergence
            declareGet(code,{full([u;d]),full(nu),full(lambda)},'getZNL__');
            declareGet(code,{full(F),full(G)},'getFG__');
        end

        lu_ww=lu(WW,[cmexfunction,'_WW.subscripts'],[cmexfunction,'_WW.values']);
        if atomicFactorization
            lu_ww=declareAlias(code,lu_ww,'lu_ww',true);
        end
        
        if ~skipAffine
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% affine scaling direction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('(%.2f sec)\n    affine direction...',etime(clock(),t2));
            t2=clock();

            dx_a=lu_ww\b_a;
            dx_a=declareAlias(code,dx_a,'dx_a__');

            fprintf('(%.2f sec)\n    updated vectores...',etime(clock(),t2));
            t2=clock();

            dU_a=dx_a(1:nU);
            dD_a=dx_a(nU+1:nU+nD);
            dZ_a=dx_a(1:nZ);
            if nF>0
                dLambda_a=dx_a(nZ+nG+nH+1:nZ+nG+nH+nF);
            else
                dLambda_a=Tzeros(nF);
            end
            
            newU_a=u+alphaPrimal*dU_a;
            newD_a=d+alphaPrimal*dD_a;
            if nX>0
                dX_a=dx_a(nU+nD+1:nZ);
                newX_a=x+alphaPrimal*dX_a;
            end
            
            if nF>0
                newLambda_a=lambda+alphaDualIneq*dLambda_a;
                
                newF_a=substitute(F,u,newU_a);
                newF_a=substitute(newF_a,d,newD_a);
                if nX>0
                    newF_a=substitute(newF_a,x,newX_a);
                end
                
                maxAlphaPrimal_a=clp(F,F_z*dZ_a);
                maxAlphaDual_a=clp(lambda,dLambda_a);
                
                % Mehrotra correction
                Mehrotra=(F_z*dZ_a).*dLambda_a./lambda;
                b_s=[-Lf_u;
                     -Lg_d;
                     -Lf_x;
                     -Lg_x;
                     -G;
                     -F-Mehrotra+muOnes./lambda];
                
                rho=(newF_a*newLambda_a)./gap;
                
                declareGet(code,{maxAlphaPrimal_a,maxAlphaDual_a},'getMaxAlphas_a__');
                declareGet(code,min(newF_a,1),'getMinF_a__');
                declareGet(code,rho,'getRho__');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% combined search direction  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('(%.2f sec)\n    combined direction...',etime(clock(),t2));
        t2=clock();
        
        dx_s=lu_ww\b_s;
        dx_s=declareAlias(code,dx_s,'dx_s__');

        dU_s=dx_s(1:nU);
        dD_s=dx_s(nU+1:nU+nD);
        dZ_s=dx_s(1:nZ);
        if nX>0
            dX_s=dx_s(nU+nD+1:nZ);
        end        
        dNu_s=dx_s(nZ+1:nZ+nG+nH);
        dNuU_s=dx_s(nZ+1:nZ+nGu);
        dNuD_s=dx_s(nZ+nGu+nH+1:nZ+nGu+nH+nGd);
        dNuUx_s=dx_s(nZ+nGu+1:nZ+nGu+nH);
        dNuDx_s=dx_s(nZ+nGu+nH+nGd+1:nZ+nGu+nH+nGd+nH);

        if nF>0
            dLambda_s=dx_s(nZ+nG+nH+1:nZ+nG+nH+nF);
        else
            dLambda_s=Tzeros(nF);
        end
        
        fprintf('(%.2f sec)\n    new vectores...',etime(clock(),t2));
        t2=clock();
        
        newU_s=u+alphaPrimal*dU_s;
        newD_s=d+alphaPrimal*dD_s;
        if nX>0
            newX_s=x+alphaPrimal*dX_s;
        end
        
        if nF>0
            maxAlphaPrimal_s=clp(F,F_z*dZ_s);
            maxAlphaDual_s=clp(lambda,dLambda_s);
            newLambda_s=lambda+alphaDualIneq*dLambda_s;
            newF_s=substitute(F,u,newU_s);
            newF_s=substitute(newF_s,d,newD_s);
            if nX>0
                newF_s=substitute(newF_s,x,newX_s);
            end
            declareGet(code,{maxAlphaPrimal_s,maxAlphaDual_s},'getMaxAlphas_s__');
            declareGet(code,min(newF_s,1),'getMinF_s__');
            if debugConvergence
                declareGet(code,{newF_s,newLambda_s},'getFLambda_s__');
            end
        end
        
    end % smallerNewtonMatrix

    if debugConvergence
        fprintf('(%.2f sec)\n    getGrad()...',etime(clock(),t2));
        t2=clock();
        
        declareGet(code,full([gradient(f,u);
                            gradient(-P1lambda*Fu,u);
                            gradient(+P1nu*Gu,u);
                            gradient(f-P1lambda*Fu+P1nu*Gu,u)]),'getGrad__');
        
        % not latent stuff
        % WWW=[Lf_uz,Lf_un,Lf_ul;
        %      Lg_dz,Lg_dn,Lg_dl;
        %      G_z,Tzeros([nG,nG]),Tzeros([nG,nF]);
        %      -F_z,Tzeros([nF,nG]),-diag(F./lambda)];
        % bb_s=-[Lf_u;
        %        Lg_d;
        %        G;
        %        -F+muOnes./lambda];
        % declareGet(code,full(WWW),'getWWW__');
        % declareGet(code,full(bb_s),'getBb_s__');
        declareGet(code,full(dZ_s),'getDz_s__');
        declareGet(code,full(dNu_s),'getDnu_s__');
        declareGet(code,full(dLambda_s),'getDlambda_s__');
    end
    
    fprintf('(%.2f sec)\n    updatePrimalDual()...',etime(clock(),t2));
    t2=clock();
    dst={u,d};
    src={newU_s,newD_s};
    if nX>0
        dst{end+1}=x;
        src{end+1}=newX_s;
    end
    if nGu>0
        dst{end+1}=P1nu;
        src{end+1}=P1nu+alphaDualEq*dNuU_s;
    end
    if nGd>0
        dst{end+1}=P2nu;
        src{end+1}=P2nu+alphaDualEq*dNuD_s;
    end
    if nH>0
        dst{end+1}=P1xnu;
        src{end+1}=P1xnu+alphaDualEq*dNuUx_s;
        dst{end+1}=P2xnu;
        src{end+1}=P2xnu+alphaDualEq*dNuDx_s;
    end
    if nFu>0
        dst{end+1}=P1lambda;
        src{end+1}=newLambda_s(1:nFu);
    end
    if nFd>0
        dst{end+1}=P2lambda;
        src{end+1}=newLambda_s(nFu+1:nF);
    end

    declareCopy(code,dst,src,'updatePrimalDual__');
    
    % declareSave after lu's to make sure previous "typical" values
    % are used by lu, prior to being overwritten by declareSave
    if allowSave
        declareSave(code,WW,'saveWW__',[cmexfunction,'_WW.subscripts']);
        %% for debug
        declareSave(code,b_s,'saveb_s__',[cmexfunction,'_b_s.subscripts']);
        if smallerNewtonMatrix
            declareSave(code,dZNu_s,'savedx_s__',[cmexfunction,'_dx_s.subscripts']);
        else
            declareSave(code,dx_s,'savedx_s__',[cmexfunction,'_dx_s.subscripts']);
        end
    end

    fprintf('(%.2f sec)\n    ',etime(clock(),t2));
    fprintf(' done ipmPDeqlat_CS symbolic computations (%.3f sec)\n',etime(clock(),t1));
    
    %profile off
    %profile viewer
    
