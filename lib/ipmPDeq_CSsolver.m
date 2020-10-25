function [varargout]=ipmPDeq_CSsolver(obj,mu0,maxIter,saveIter)
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

    FUNCTION__='ipmPDeq_CSsolver';

    if nargin<2
        mu0=1;
    end
    if nargin<3
        maxIter=200;
    end
    if nargin<4
        saveIter=-1;
    end

    function printf2(varargin)
        if obj.verboseLevel>=2
            fprintf(varargin{:});
        end
    end

    function printf3(varargin)
        if obj.verboseLevel>=3
            fprintf(varargin{:});
        end
    end

    iter=0;

    initPrimal__(obj);

    if obj.scaleCost
        scaleCost__(obj);
        desiredDualityGap=getScale4Cost__(obj)*obj.desiredDualityGap;
    else
        desiredDualityGap=obj.desiredDualityGap;
    end

    if obj.nF>0
        if obj.scaleInequalities
            scaleIneq__(obj);
        end
        alphaPrimal=0;alphaDualEq=0;alphaDualIneq=0;
        mu=mu0;
        setMu__(obj,mu);
        muMin=desiredDualityGap/obj.nF/2;
    end

    printf2('%s.m (coupledAlphas=%d,skipAffine=%d,delta=%g): %d primal variables (%d+%d+%d), %d eq. constr., %d ineq. constr.\n',...
            FUNCTION__,obj.coupledAlphas,obj.skipAffine,obj.delta,obj.nZ,obj.nU,obj.nD,obj.nX,obj.nG,obj.nF);
    if obj.verboseLevel>=3
        headers='Iter   cost1      cost2      |grad|     |eq|    inequal     dual      gap       mu      alphaA     sigma   alphaP    alphaDI     alphaDE  time [ms]\n';
        fprintf(headers);
    end
    if obj.nF>0
        printf3('%3d: <-maxIter       tol->%10.2e%10.2e                    %10.2e%10.2e\n',...
                maxIter,obj.gradTolerance,obj.equalTolerance,desiredDualityGap,muMin);
    else
        printf3('%3d: <-maxIter       tol->%10.2e%10.2e\n',maxIter,obj.gradTolerance,obj.equalTolerance);
    end

    dt0=clock();

    if obj.nF>0
        initDualIneq__(obj);
        if obj.debugConvergence
            [F_,l_]=getFLambda__(obj);
            k=find(F_<1/sqrt(obj.debugConvergenceThreshold));
            if ~isempty(k)
                printf2('%3d: ATTENTION: initial ineq < %10.2e for %4d entries - scale optimization variables or add constraint\n',...
                        iter,1/sqrt(obj.debugConvergenceThreshold),length(k));
                if obj.verboseLevel>=4
                    for ii=k(:)'
                        fprintf('\t ineq %4d: %9.2e\n',ii,full(F_(ii)));
                    end
                end
            end
        end
    end

    if obj.nG>0
        initDualEq__(obj);
    end

    while (1)
        iter=iter+1;
        if obj.verboseLevel>=3
            if mod(iter,50)==0
                fprintf(headers);
            end
            fprintf('%3d:',iter);
            dt1=clock();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Check exit conditions %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        if iter > maxIter
            printf3('maximum # iterations (%d) reached.\n',maxIter);
            status = 8;
            break;
        end

        if obj.verboseLevel>=3 || obj.debugConvergence
            [f,g]=getfg__(obj);
        end
        if obj.verboseLevel>=3
            fprintf('%11.3e%11.3e',full(f),full(g));
        end

        norminf_grad=getNorminf_Grad__(obj);
        printf3('%10.2e',full(norminf_grad));

        if obj.debugConvergence
            if norminf_grad>obj.debugConvergenceThreshold
                grad_=getGrad__(obj);
                k=find(abs(grad_)>obj.debugConvergenceThreshold);
                if ~isempty(k)
                    printf2('\n%3d: ATTENTION: |grad| > %10.2e for %4d entries - cost is probably too large or primal variables too small\n',...
                            iter,obj.debugConvergenceThreshold,length(k));
                    for ii=k(:)'
                        printf2('\t grad (primal var %3d): %9.2e\n',ii,grad_(ii));
                    end
                end
            end
        end

        if isnan(norminf_grad)
            printf2('  -> failed to invert hessian\n');
            status = 4;
            break;
        end

        if obj.nG>0
            norminf_eq=getNorminf_G__(obj);
            printf3('%10.2e',full(norminf_eq));
        else
            printf3('    -eq-  ');
        end

        if obj.nF>0
            [gap,ineq,dual]=getGapMinFMinLambda__(obj);
            printf3('%10.2e%10.2e%10.2e',full(ineq),full(dual),full(gap));
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
            printf3('   -ineq-    -dual-    -gap-  ');
        end

        if norminf_grad<=obj.gradTolerance && ...
                (obj.nF==0 || gap<=desiredDualityGap) && ...
                (obj.nG==0 || norminf_eq<=obj.equalTolerance)
            printf2('  -> clean exit\n');
            status = 0;
            break;
        end

        if obj.nF>0
            printf3('%10.2e',mu);
        else
            printf3('   -mu-   ');
        end

        if obj.debugConvergence
            [Lfu,Lgd]=getLf1__(obj);
            [Lfuu,Lgdd]=getLf2__(obj);
            euu=eig(Lfuu);
            tol=1e-9;
            mpu=sum(euu>tol);
            mnu=sum(euu<-tol);
            edd=eig(Lgdd);
            mpd=sum(edd>tol);
            mnd=sum(edd<-tol);
            fprintf('\n   |Lfu|=%g, |Lgd|=%g, |Lfuu-Lfuu''|=%g, |Lgdd-Lgdd''|=%g, eig(Lfuu) in [%g,%g] (%d >0,%d =0,%d <0), eig(Lgdd) in [%g,%g] (%d >0,%d =0,%d <0)\n',...
                    norm(Lfu,'fro'),norm(Lgd,'fro'),...
                    norm(Lfuu-Lfuu','fro'),norm(Lgdd-Lgdd','fro'),...
                    min(euu),max(euu),mpu,length(euu)-mpu-mnu,mnu,min(edd),max(edd),mpd,length(edd)-mpd-mnd,mnd);
        end


        if obj.nF==0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  NO INEQUALITY CONSTRAINTS %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            setAlphaPrimal__(obj,obj.alphaMax);
            if obj.nG>0
                setAlphaDualEq__(obj,obj.alphaMax);
            end
            printf3('  -alphaA-  -sigma- ');
            printf3('%10.2e                   ',obj.alphaMax);

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
                printf3(' -alphaA-  -sigma-');
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Affine search direction %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                [alphaPrimal,alphaDualIneq]=getMaxAlphas_a__(obj);

                alphaMax = min([obj.alphaMax,alphaPrimal,alphaDualIneq]);

                if (alphaMax >= obj.alphaMin)
                    % try max
                    alphaPrimal=alphaMax;
                    setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_a__(obj);
                    if (ineq<0)
                        % try min
                        alphaPrimal=obj.alphaMin;
                        setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_a__(obj);
                        if (ineq>0)
                            % try between min and max
                            alphaPrimal = alphaMax*.95;
                            while alphaPrimal >= obj.alphaMin
                                setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_a__(obj);
                                if (ineq>=0)
                                    break;
                                end
                                alphaPrimal=alphaPrimal/2;
                            end
                            if (alphaPrimal < obj.alphaMin)
                                alphaPrimal = 0;
                                setAlphaPrimal__(obj,alphaPrimal);
                            end
                        else
                            alphaPrimal = 0;
                            setAlphaPrimal__(obj,alphaPrimal);
                        end
                    end
                else
                    alphaPrimal = 0;
                    setAlphaPrimal__(obj,alphaPrimal);
                end
                setAlphaDualIneq__(obj,alphaPrimal);
                printf3('%10.2e',full(alphaPrimal));

                % update mu based on sigma, but this only seems to be safe for:
                % 1) 'long' newton steps in the affine direction
                % 2) equality constraints fairly well satisfied (perhaps not very important)
                % 3) small gradient
                %th_grad=norminf_grad<=max(1e-1,1e2*obj.gradTolerance);
                th_eq=(obj.nG==0) || norminf_eq<=1e-3 || norminf_eq<=1e2*obj.equalTolerance;
                if alphaPrimal>obj.alphaMax/2 && th_eq %&& th_grad
                    sigma=full(getRho__(obj));
                    if (sigma>1) sigma=1; end
                    if (sigma<0) sigma=0; end
                    if obj.delta==2
                        sigma=sigma*sigma;
                    else
                        sigma=sigma*sigma*sigma;
                    end
                    printf3('%10.2e',sigma);
                    mu = full(max(sigma*gap/obj.nF,muMin));
                    setMu__(obj,mu);
                else
                    printf3('  -sigma- ');
                end
            end  % obj.skipAffine==1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Combined search direction %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if obj.debugConvergence
                setAlphaPrimal__(obj,1);
                if obj.nG>0
                    setAlphaDualEq__(obj,1);
                end
                if obj.nF>0
                    setAlphaDualIneq__(obj,1);
                end
                [newF_s,newLambda_s]=getFLambda_s__(obj);
                newmu = mu;
            end

            [alphaPrimal,alphaDualIneq]=getMaxAlphas_s__(obj);

            if obj.coupledAlphas && alphaDualIneq<alphaPrimal
                alphaPrimal=alphaDualIneq;
            end

            alphaPrimal = .99 * alphaPrimal;

            alphaMax = min(alphaPrimal,obj.alphaMax);

            if (alphaMax >= obj.alphaMin)
                % try max
                alphaPrimal=alphaMax/.99;
                setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_s__(obj);
                if isnan(ineq)
                    printf2('  -> failed to invert hessian\n');
                    status = 4;
                    break;
                end
                if (ineq>0)
                    % recheck just to be safe in case not convex
                    alphaPrimal = .99 * alphaPrimal;
                    setAlphaPrimal__(obj,alphaPrimal);ineq1=getMinF_s__(obj);
                end
                if ineq<=0 || ineq1<=ineq/10
                    % try min
                    alphaPrimal=obj.alphaMin/.99;
                    setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_s__(obj);
                    if (ineq>0)
                        % try between min and max
                        alphaPrimal=alphaMax*.95;
                        while alphaPrimal >= obj.alphaMin
                            setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_s__(obj);
                            if (ineq>0)
                                % backtrace just a little
                                alphaPrimal = .99 * alphaPrimal;
                                % recheck just to be safe in case not convex
                                setAlphaPrimal__(obj,alphaPrimal);ineq1=getMinF_s__(obj);
                                if (ineq1>ineq/10)
                                    break;
                                end
                            end
                            alphaPrimal=alphaPrimal/2;
                        end
                        if (alphaPrimal < obj.alphaMin)
                            alphaPrimal = 0;
                            setAlphaPrimal__(obj,alphaPrimal);
                        end
                    else
                        alphaPrimal = 0;
                        setAlphaPrimal__(obj,alphaPrimal);
                    end
                end
            else
                alphaPrimal = 0;
                setAlphaPrimal__(obj,alphaPrimal);
            end

            if obj.coupledAlphas
                alphaDualEq=alphaPrimal;
                alphaDualIneq=alphaPrimal;
            else
                alphaDualIneq = .99 * alphaDualIneq;
                if alphaDualIneq>obj.alphaMax
                    alphaDualIneq = obj.alphaMax;
                end
                alphaDualEq = obj.alphaMax;
            end

            if obj.nG>0
                setAlphaDualEq__(obj,alphaDualEq);
            end
            setAlphaDualIneq__(obj,alphaDualIneq);
            updatePrimalDual__(obj);

            if obj.nG>0
                printf3('%10.2e %10.2e %10.2e',full(alphaPrimal),full(alphaDualIneq),full(alphaDualEq));
            else
                printf3('%10.2e %10.2e   -eq-    ',full(alphaPrimal),full(alphaDualIneq));
            end

            if obj.skipAffine==1
                % More aggressive if
                % 1) 'long' newton steps in the affine direction
                % 2) small gradient
                % 3) equality constraints fairly well satisfied
                % (2+3 mean close to the central path)
                %th_grad=norminf_grad<=max(1e-1,1e2*obj.gradTolerance);
                %th_eq=(obj.nG==0) || (norminf_eq<=max(1e-3,1e2*obj.equalTolerance));
                th_grad=            norminf_grad<=max(1e-3,1e0*obj.gradTolerance);
                th_eq=(obj.nG==0) || (norminf_eq<=max(1e-5,1e0*obj.equalTolerance));
                if alphaPrimal>obj.alphaMax/2 && th_grad && th_eq
                    %mu = max(mu*obj.muFactorAggressive,muMin);
                    mu = max(muMin,min(obj.muFactorAggressive*mu,mu^1.5));
                    setMu__(obj,mu);
                    printf3(' * ');
                else
                    if alphaPrimal<.1
                        mu=min(mu0,1.1*mu);
                        setMu__(obj,mu);
                        initDualIneq__(obj);
                        printf3('^');
                    else
                        mu=max(mu*obj.muFactorConservative,muMin);
                        setMu__(obj,mu);
                        printf3('v');
                    end
                    if th_grad
                        printf3('g');
                    else
                        printf3(' ');
                    end
                    if th_eq
                        printf3('e');
                    else
                        printf3(' ');
                    end
                end
            end

            % if no motion, slowly increase mu
            if (alphaPrimal<obj.alphaMin && alphaDualIneq<obj.alphaMin && alphaDualEq<obj.alphaMin)
                mu = max(mu/obj.muFactorConservative,muMin);
                setMu__(obj,mu);
            end

        end  % if obj.nF==0

        if obj.verboseLevel>=3
            dt1=etime(clock(),dt1);
            fprintf('%8.1fms\n',dt1*1e3);
        end

        %%%%%%%%%%%%%%%%%%%%%%%
        %% Debug small alpha %%
        %%%%%%%%%%%%%%%%%%%%%%%

        if obj.debugConvergence
            if obj.nF>0 && alphaPrimal<obj.alphaMax/5 && isfinite(obj.debugConvergenceThreshold)
                kk=find(newF_s<=0 | newLambda_s<=0);
                fprintf('%3d: ATTENTION: alphaPrimal = %8.2e < %8.2e due to %4d entries (lambda * ineq)\n',iter,full(alphaPrimal),obj.alphaMax/5,length(kk));
                if obj.verboseLevel>=4
                    for ii=kk(:)'
                        fprintf('\t ineq %4d: %8.2e * %8.2e = %8.2e (%8.2e) -> %9.2e * %9.2e = %9.2e (%8.2e)',...
                                ii,full(oldLambda(ii)),full(oldF(ii)),full(oldLambda(ii)*oldF(ii)),oldmu,...
                                full(newLambda_s(ii)),full(newF_s(ii)),full(newLambda_s(ii)*newF_s(ii)),newmu);
                        if (newF_s(ii)<0)
                            fprintf(' -- probably inequality too large (needs scaling)\n');
                        else
                            fprintf(' -- probably inequality too large or mu too small\n');
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%
        %% Check scaling %%
        %%%%%%%%%%%%%%%%%%%

        if obj.debugConvergence
            %% Check scaling for optimization variables
            [u_,d_]=getUD__(obj);
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

            k=find(abs(d_)>obj.debugConvergenceThreshold);
            if ~isempty(k)
                fprintf('%3d: ATTENTION: abs(d) > %10.2e > %10.2e for %4d entries - scale optimization variables or add constraint\n',...
                        iter,obj.debugConvergenceThreshold,length(k));
                for ii=k(:)'
                    fprintf('\t primal var %4d: %9.2e\n',ii,d_(ii));
                end
            end
            if max(abs(d_))<1/obj.debugConvergenceThreshold
                %disp(d_')
                fprintf('%3d: ATTENTION: all abs(d) < %10.2e - scale optimization variables\n',...
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
                        fprintf('\t ineq %4d: %9.2e\n',ii,full(F_(ii)));
                    end
                end
                if max(abs(l_))<1/obj.debugConvergenceThreshold
                    %disp(l_')
                    fprintf('%3d: ATTENTION: all abs(lambda) < %10.2e (%10.2e<=abs(F)<=%10.2e) - either all inequalities inactive or scale inequality constraint\n',...
                            iter,1/obj.debugConvergenceThreshold,min(abs(F_)),max(abs(F_)));
                    if obj.verboseLevel>=4
                        disp([l_,F_])
                    end
                end
            end
        end

    end % while(1)

    if obj.debugConvergence
        [Lfu,Lgd]=getLf1__(obj);
        [Lfuu,Lgdd]=getLf2__(obj);
        euu=eig(Lfuu);
        tol=1e-9;
        mpu=sum(euu>tol);
        mnu=sum(euu<-tol);
        edd=eig(Lgdd);
        mpd=sum(edd>tol);
        mnd=sum(edd<-tol);
        fprintf('   |Lfu|=%g, |Lgd|=%g, |Lfuu-Lfuu''|=%g, |Lgdd-Lgdd''|=%g, eig(Lfuu) in [%g,%g] (%d >0,%d =0,%d <0), eig(Lgdd) in [%g,%g] (%d >0,%d =0,%d <0)\n',...
                norm(Lfu,'fro'),norm(Lgd,'fro'),...
                norm(Lfuu-Lfuu','fro'),norm(Lgdd-Lgdd','fro'),...
                min(euu),max(euu),mpu,length(euu)-mpu-mnu,mnu,min(edd),max(edd),mpd,length(edd)-mpd-mnd,mnd);
    end

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
            if (gap>desiredDualityGap)
                status=bitor(status,64);
            end
            if (mu>muMin)
                status=bitor(status,128);
            end
            if (alphaPrimal<=obj.alphaMin && alphaDualIneq<=obj.alphaMin && alphaDualEq<=obj.alphaMin)
                status=bitor(status,1792); % (256|512|1024);
            elseif (alphaPrimal<=.1 && alphaDualIneq<=.1 && alphaDualEq<=.1)
                status=bitor(status,1536); % (512|1024);
            elseif (alphaPrimal<=.5 && alphaDualIneq<=.5 && alphaDualEq<=.5)
                status=bitor(status,1024);
            end
        end
    end

    time=etime(clock(),dt0);
    if obj.verboseLevel>=2
        [f,g]=getfg__(obj);
        if obj.nG>0 && status<8 % when status>=8 this has been already been computed
            norminf_eq=getNorminf_G__(obj);
        end
        if obj.nF>0 && status<8 % when status>=8 this has been already been computed
            [gap,ineq,dual]=getGapMinFMinLambda__(obj);
        end

        if (status)
            fprintf('%3d:status=0x%s ',iter,dec2hex(status));
            sep='(';
            if bitand(status,16)
                fprintf("%clarge gradient",sep);
                sep=',';
            end
            if bitand(status,32)
                fprintf("%cbad equality const.",sep);
                sep=',';
            end
            if bitand(status,64)
                fprintf("%clarge duality gap",sep);
                sep=',';
            end
            if bitand(status,128)
                fprintf("%clarge mu",sep);
                sep=',';
            end
            if bitand(status,256)
                fprintf("%calpha negligible",sep);
                sep=',';
            elseif bitand(status,512)
                fprintf("%calpha<.1",sep);
                sep=',';
            elseif bitand(status,1024)
                fprintf("%calpha<.5",sep);
                sep=',';
            end
            fprintf(')\n                ');
        else
            fprintf('%3d:status=0x%s, ',iter,dec2hex(status));
        end

        fprintf('cost=%13.5e,%13.5e ',full(f),full(g));
        norminf_grad=getNorminf_Grad__(obj);
        fprintf('|grad|=%10.2e',full(norminf_grad));
        if obj.nG>0
            fprintf(', |eq|=%10.2e',full(norminf_eq));
        end
        if obj.nF>0
            fprintf(', ineq=%10.2e,\n              dual=%10.2e, gap=%10.2e, last alpha=%10.2e, last mu=%10.2e',full(ineq),full(dual),full(gap),full(alphaPrimal),mu);
        end
        fprintf(' (%.1fms,%.2fms/iter)\n',time*1e3,time/iter*1e3);
    end

    if nargout==1
        varargout{1}.status=status;
        varargout{1}.iter=iter;
        varargout{1}.time=time;
    else
        varargout={status,iter,time};
    end

end