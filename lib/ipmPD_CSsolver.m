function [status,iter,time]=ipmPD_CSsolver(obj,mu0,maxIter,saveIter)
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

    FUNCTION__='ipmPD_CSsolver';
    
    function printf1(varargin)
        if obj.verboseLevel>=2
            fprintf(varargin{:});
        end
    end
    
    function printf2(varargin)
        if obj.verboseLevel>=3
            fprintf(varargin{:});
        end
    end
    
    iter=0;
    
    if obj.nF>0
        mu=mu0;
        alpha=0;
        muMin=obj.desiredDualityGap/obj.nF/2;
    end 
    
    printf1('%s.m (skipAffine=%d,delta=%g): %d primal variable, %d equality constraints, %d inequality constraints\n',FUNCTION__,obj.skipAffine,obj.delta,obj.nU,obj.nG,obj.nF);
    printf2('Iter   cost      |grad|     |eq|     inequal     dual      gap       mu      alphaA    sigma     alphaS   time [ms]\n');
    if obj.nF>0
        printf2('%3d:<-mx tol-> %10.2e%10.2e                    %10.2e%10.2e\n',...
                maxIter,obj.gradTolerance,obj.equalTolerance,obj.desiredDualityGap,muMin);
    else
        printf2('%3d:<-mx tol-> %10.2e%10.2e\n',maxIter,obj.gradTolerance,obj.equalTolerance);
    end
    
    dt0=clock();
    
    %initPrimalDual__(obj);
    initPrimal__(obj);
    
    if obj.nF>0
        setMu__(obj,mu);
    end
    
    if obj.nG>0
        initDualEq__(obj);
    end
    
    if obj.nF>0
        initDualIneq__(obj);
        if obj.debugConvergence 
            [F_,l_]=getFLambda__(obj);
            k=find(F_<1/sqrt(obj.debugConvergenceThreshold));
            if ~isempty(k)
                printf1('%3d: ATTENTION: initial ineq < %10.2e for %4d entries - scale optimization variables or add constraint\n',...
                        iter,1/sqrt(obj.debugConvergenceThreshold),length(k));
                for ii=k(:)'
                    printf1('\t ineq %4d: %9.2e\n',ii,F_(ii));
                end
            end
        end
    end
    
    if obj.debugConvergence 
        lastJ=inf;
    end

    while (1) 
        if obj.verboseLevel>=3
            dt1=clock();
        end
        
        iter=iter+1;
        printf2('%3d:',iter);
        
        if iter > maxIter 
            printf2('maximum # iterations (%d) reached.\n',maxIter);
            status = 8;
            break; 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Check exit conditions %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if obj.verboseLevel>=3
            J=getJ__(obj);
            fprintf('%11.3e',full(J));
        end

        if obj.debugConvergence 
            if J>obj.debugConvergenceThreshold
                printf1('\n%3d: ATTENTION: cost > %10.2e - cost is probably too large\n',...
                            iter,obj.debugConvergenceThreshold);
            end
        end

        norminf_grad=getNorminf_Grad__(obj);
        printf2('%10.2e',norminf_grad);
        
        if obj.debugConvergence 
            if norminf_grad>obj.debugConvergenceThreshold
                grad_=getGrad__(obj);
                k=find(abs(grad_)>obj.debugConvergenceThreshold);
                if ~isempty(k)
                    printf1('\n%3d: ATTENTION: |grad| > %10.2e for %4d entries - cost is probably too large or primal variables too small\n',...
                            iter,obj.debugConvergenceThreshold,length(k));
                    for ii=k(:)'
                        printf1('\t grad (primal var %3d): %9.2e\n',ii,grad_(ii));
                    end
                end
            end
        end

        if obj.debugConvergence 
            g_th=max(abs(J-lastJ)*obj.debugConvergenceThreshold,10*obj.gradTolerance);
            if abs(J-lastJ)<1/obj.debugConvergenceThreshold && norminf_grad>g_th
                grad_=getGrad__(obj);
                k=find(abs(grad_)>g_th);
                if ~isempty(k)
                    printf1('\n%3d: ATTENTION: |Delta J| < %10.2e, but |grad| > %10.2e for %4d entries - poorly conditioned hessian\n',...
                            iter,1/obj.debugConvergenceThreshold,g_th,length(k));
                    for ii=k(:)'
                        printf1('\t grad %3d: %9.2e\n',ii,grad_(ii));
                    end
                end
            end
            lastJ=J;
        end
        
        if isnan(norminf_grad) 
            printf2('  -> failed to invert hessian\n');
            status = 4;
            break;
        end
        
        if obj.nG>0
            norminf_eq=getNorminf_G__(obj);
            if obj.verboseLevel>=3
                printf2('%10.2e',norminf_eq);
            end
        else
            printf2('    -eq-  ');
        end
        
        if obj.nF>0
            [gap,ineq,dual]=getGapMinFMinLambda__(obj);
            printf2('%10.2e%10.2e%10.2e',ineq,dual,gap);
            if (ineq<=0) 
                printf2('  -> (primal) variables violate constraints\n');
                status = 1;
                break;
            end
            if (dual<=0) 
                printf2('  -> negative value for dual variables\n');
                    status = 2;
                    break;
            end
        else
            printf2('   -ineq-    -dual-    -gap-  ');
        end
        
        if norminf_grad<=obj.gradTolerance && ...
                (obj.nF==0 || gap<=obj.desiredDualityGap) && ...
                (obj.nG==0 || norminf_eq<=obj.equalTolerance)
            printf2('  -> clean exit\n');
            status = 0;
            break;
        end
        
        if obj.nF>0
            printf2('%10.2e',mu);
        else
            printf2('   -mu-   ');
        end
        
        % WW__=getWW__(obj)
        % b_s__=getb_s__(obj)
        % dx_s__=getDx_s__(obj)
        % u__=getU__(obj)
        % lambda__=getLambda__(obj)
        
        if obj.nF==0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  NO INEQUALITY CONSTRAINTS %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            setAlpha__(obj,obj.alphaMax);
            printf2('  -alphaA-  -sigma- ');
            printf2('%10.2e',obj.alphaMax);
            
            updatePrimalDual__(obj);
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  WITH INEQUALITY CONSTRAINTS %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.debugConvergence
                [oldF,oldLambda]=getFLambda__(obj);
                oldmu=mu;
            end
            if obj.skipAffine==1
                printf2(' -alphaA-  -sigma-');
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Affine search direction %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [primalAlpha,dualAlpha]=getAlphas_a__(obj);
                
                alpha = min(primalAlpha,dualAlpha);
                alphaMax = min(alpha,obj.alphaMax);
                %fprintf('\nAlphaPrimal_a = %10.3e, AlphaDual_a = %10.3e, alpha = %10.3e, alphaMax = %10.3e, alphaMin = %10.3e\n',primalAlpha,dualAlpha,alpha,alphaMax,obj.alphaMin);
                
                if (alphaMax >= obj.alphaMin) 
                    % try max
                    alpha=alphaMax;
                    setAlpha__(obj,alpha);ineq=getMinF_a__(obj);
                    if (ineq<0) 
                        % try min
                        alpha=obj.alphaMin;
                        setAlpha__(obj,alpha);ineq=getMinF_a__(obj);
                        if (ineq>0) 
                            % try between min and max
                            alpha = alphaMax*.95;
                            while alpha >= obj.alphaMin
                                setAlpha__(obj,alpha);ineq=getMinF_a__(obj);
                                if (ineq>=0) 
                                    break; 
                                end
                                alpha=alpha/2;
                            end
                            if (alpha < obj.alphaMin) 
                                alpha = 0;
                                setAlpha__(obj,alpha);
                            end
                        else 
                            alpha = 0;
                            setAlpha__(obj,alpha);
                        end
                    end
                else 
                    alpha = 0;
                    setAlpha__(obj,alpha);
                end
                printf2('%10.2e',alpha);
                
                % update mu based on sigma, but this only seems to be safe for:
                % 1) 'long' newton steps in the affine direction
                % 2) equality constraints fairly well satisfied (perhaps not very important)
                % 3) small gradient
                %th_grad=norminf_grad<=max(1e-1,1e2*obj.gradTolerance);
                th_eq=(obj.nG==0) || (norminf_eq<=max(1e-3,1e2*obj.equalTolerance));
                if alpha>obj.alphaMax/2 && th_eq %&& th_grad 
                    sigma=getRho__(obj);
                    if (sigma>1) sigma=1; end
                    if (sigma<0) sigma=0; end
                    if obj.delta==2
                        sigma=sigma*sigma;
                    else
                        sigma=sigma*sigma*sigma;
                    end
                    printf2('%10.2e',sigma);
                    mu=max(sigma*gap/obj.nF,muMin);
                    setMu__(obj,mu); 
                else 
                    printf2('  -sigma- ');
                end
            end  % obj.skipAffine==1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Combined search direction %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.debugConvergence
                setAlpha__(obj,1);
                [newF_s,newLambda_s]=getFLambda_s__(obj);
                newmu=mu;
            end

            [primalAlpha,dualAlpha]=getAlphas_s__(obj);

            alpha = .99*min(primalAlpha,dualAlpha);
            alphaMax = min(alpha,obj.alphaMax);
            %fprintf('\nAlphaPrimal_s = %10.3e, AlphaDual_s = %10.3e, alpha = %10.3e, alphaMax = %10.3e, alphaMin = %10.3e\n',primalAlpha,dualAlpha,alpha,alphaMax,obj.alphaMin);
            
            if (alphaMax >= obj.alphaMin) 
                % try max
                alpha=alphaMax/.99;
                setAlpha__(obj,alpha);ineq=getMinF_s__(obj);
                %printf(' minF(maxAlpha)=%10.3e ',ineq);
                if isnan(ineq) 
                    printf2('  -> failed to invert hessian\n');
                    status = 4;
                    break;
                end
                if (ineq>0) 
                    alpha = alpha*.99;
                    setAlpha__(obj,alpha);ineq1=getMinF_s__(obj);
                    if ineq1>ineq/10
                        updatePrimalDual__(obj);
                    end
                end
                if ineq<=0 || ineq1<=ineq/10
                    % try min
                    alpha=obj.alphaMin/.99;
                    setAlpha__(obj,alpha);ineq=getMinF_s__(obj);
                    %printf(' minF(minAlpha)=%10.3e ',ineq);
                    if (ineq>0) 
                        % try between min and max
                        alpha=alphaMax*.95;
                        while alpha >= obj.alphaMin
                            setAlpha__(obj,alpha);ineq=getMinF_s__(obj);
                            %printf(' minF(%g)=%10.3e ',alpha,ineq);
                            if (ineq>0) 
                                % backtrace just a little
                                alpha = alpha*.99; 
                                % recheck just to be safe in case not convex
                                setAlpha__(obj,alpha);ineq1=getMinF_s__(obj);
                                if (ineq1>ineq/10)
                                    updatePrimalDual__(obj);
                                    break; 
                                end
                            end
                            alpha=alpha/2;
                        end
                        if (alpha < obj.alphaMin) 
                            alpha = 0;
                            setAlpha__(obj,alpha);
                        end
                    else 
                        alpha = 0;
                        setAlpha__(obj,alpha);
                    end
                end
            else 
                alpha = 0;
                setAlpha__(obj,alpha);
            end
            
            printf2('%10.2e',alpha);

            if obj.skipAffine==1
                % More aggressive if
                % 1) 'long' newton steps in the affine direction
                % 2) small gradient
                % 3) equality constraints fairly well satisfied
                % (2+3 mean close to the central path)
                th_grad=norminf_grad<=max(1e-1,1e2*obj.gradTolerance);
                th_eq=(obj.nG==0) || (norminf_eq<=max(1e-3,1e2*obj.equalTolerance));
                if alpha>obj.alphaMax/2 && th_grad && th_eq
                    mu = max(mu*obj.muFactorAggressive,muMin);
                    setMu__(obj,mu); 
                    printf2(' * ');
                else 
                    if alpha<.1
                        mu=min(1e2,1.25*mu);
                        setMu__(obj,mu); 
                        initDualIneq__(obj);
                        printf2('^');
                    else
                        mu=max(mu*obj.muFactorConservative,muMin);
                        setMu__(obj,mu); 
                        printf2('v');
                    end
                    if th_grad
                        printf2('g');
                    else 
                        printf2(' ');
                    end
                    if th_eq
                        printf2('e');
                    else
                        printf2(' ');
                    end
                end
            end
            
            % if no motion, slowly increase mu
            if (alpha<obj.alphaMin) 
                mu=max(mu/obj.muFactorConservative,muMin);
                setMu__(obj,mu); 
            end
            
        end
        if obj.verboseLevel>=3
            dt1=etime(clock(),dt1);
            fprintf('%8.1fms\n',dt1*1e3);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %% Debug small alpha %%
        %%%%%%%%%%%%%%%%%%%%%%%
        
        if obj.debugConvergence
            if obj.nF>0 && alpha<obj.alphaMax/5
                kk=find(newF_s<=0 | newLambda_s<=0);
                fprintf('%3d: ATTENTION: alpha = %8.2e < %8.2e due to %4d entries (lambda * ineq)\n',iter,alpha,obj.alphaMax/5,length(kk));
                for ii=kk(:)'
                    fprintf('\t ineq %4d: %8.2e * %8.2e = %8.2e (%8.2e) -> %9.2e * %9.2e = %9.2e (%8.2e)',...
                            ii,oldLambda(ii),oldF(ii),oldLambda(ii)*oldF(ii),oldmu,...
                            newLambda_s(ii),newF_s(ii),newLambda_s(ii)*newF_s(ii),newmu);
                    if (newF_s(ii)<0)
                        fprintf(' -- probably inequality too large (needs scaling)\n');
                    else
                        fprintf(' -- probably inequality too large or mu too small\n');
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%
        %% Check scaling %%
        %%%%%%%%%%%%%%%%%%%
        
        if obj.debugConvergence
            %% Check scaling for optimization variables
            u_=getU__(obj);
            k=find(abs(u_)>obj.debugConvergenceThreshold);
            if ~isempty(k)
                fprintf('%3d: ATTENTION: abs(u) > %10.2e for %4d entries - scale optimization variables or add constraint\n',...
                        iter,obj.debugConvergenceThreshold,length(k));
                for ii=k(:)'
                    fprintf('\t primal var %4d: %9.2e\n',ii,u_(ii));
                end
            end
            if max(abs(u_))<1/obj.debugConvergenceThreshold
                %disp(u_')
                fprintf('%3d: ATTENTION: all abs(u) < %10.2e - scale optimization variables\n',...
                        iter,1/obj.debugConvergenceThreshold);
            end
            
            if obj.nG>0
                %% Check scaling for equality constraints
                [G_,nu_]=getGNu__(obj);
                k=find(abs(nu_)>obj.debugConvergenceThreshold);
                if ~isempty(k)
                    fprintf('%3d: ATTENTION: abs(nu) > %10.2e for %4d entries - scale equality constraints\n',...
                            iter,obj.debugConvergenceThreshold,length(k));
                    for ii=k(:)'
                        fprintf('\t eq %4d: %9.2e\n',ii,nu_(ii));
                    end
                end
                if max(abs(nu_))<1/obj.debugConvergenceThreshold
                    %disp(nu_')
                    fprintf('%3d: ATTENTION: all abs(nu) < %10.2e - scale equality constraints\n',...
                            iter,1/obj.debugConvergenceThreshold);
                end
            end

            if obj.nF>0
                %% Check scaling for inequalities constraints
                [F_,l_]=getFLambda__(obj);
                k=find(abs(l_)>obj.debugConvergenceThreshold);
                if ~isempty(k)
                    fprintf('%3d: ATTENTION: abs(lambda) > %10.2e for %4d entries - scale inequality constraint\n',...
                            iter,obj.debugConvergenceThreshold,length(k));
                    for ii=k(:)'
                        fprintf('\t ineq %4d: %9.2e\n',ii,l_(ii));
                    end
                end
                k=find(abs(F_)>obj.debugConvergenceThreshold);
                if ~isempty(k)
                    fprintf('%3d: ATTENTION: ineq > %10.2e for %4d entries - scale inequality constraint\n',...
                            iter,obj.debugConvergenceThreshold,length(k));
                    for ii=k(:)'
                        fprintf('\t ineq %4d: %9.2e\n',ii,F_(ii));
                    end
                end
                if max(abs(l_))<1/obj.debugConvergenceThreshold
                    %disp(l_')
                    fprintf('%3d: ATTENTION: all abs(lambda) < %10.2e (%10.2e<=abs(F)<=%10.2e) - either all inequalities inactive or scale inequality constraint\n',...
                            iter,1/obj.debugConvergenceThreshold,min(abs(F_)),max(abs(F_)));
                    disp([l_,F_])
                end
            end
        end        

    end % while(1)
    
    if status == 8
        norminf_grad=getNorminf_Grad__(obj);
        if (norminf_grad>obj.gradTolerance) 
            status=bitor(status,16);
        end
        if (obj.nG>0) 
            norminf_eq=getNorminf_G__(obj);
            if (norminf_eq>obj.equalTolerance)
                status=bitor(status,32);
            end
        end
        if (obj.nF>0)
            [gap,ineq,dual]=getGapMinFMinLambda__(obj);
            if (gap>obj.desiredDualityGap)
                status=bitor(status,64);
            end
            if (mu>muMin)
                status=bitor(status,128);
            end
            if (alpha<=obj.alphaMin)
                status=bitor(status,1792); % (256|512|1024);
            elseif (alpha<=.1)
                status=bitor(status,1536); % (512|1024);
            elseif (alpha<=.5)
                status=bitor(status,1024);
            end
        end
    end
    
    time=etime(clock(),dt0);
    if obj.verboseLevel>=2
        J=getJ__(obj);
        if obj.nG>0 && status<8 % when status>=8 this has been already been computed
            norminf_eq=getNorminf_G__(obj);
        end
        if obj.nF>0 && status<8 % when status>=8 this has been already been computed
            [gap,ineq,dual]=getGapMinFMinLambda__(obj);
        end
        
        printf1('%3d:status=0x%s, ',iter,dec2hex(status));
        printf1('cost=%13.5e, ',full(J));
        norminf_grad=getNorminf_Grad__(obj);
        printf1('|grad|=%10.2e',norminf_grad);
        if obj.nG>0
            printf1(', |eq|=%10.2e',norminf_eq);
        end
        if obj.nF>0
            printf1(', ineq=%10.2e,\n              dual=%10.2e, gap=%10.2e, last alpha=%10.2e',ineq,dual,gap,alpha);
        end
        printf1(' (%.1fms,%.2fms/iter)\n',time*1e3,time/iter*1e3);
    end
    
end