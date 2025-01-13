function varargout=ipmPD_CSsolver(obj,mu0,maxIter,saveIter,addEye2Hessian)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    FUNCTION__='ipmPD_CSsolver';

    %alphas=[];

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize addEye2Hessian %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if obj.setAddEye2Hessian
        addEye2HessianUMAX=1e2;
        addEye2HessianEqMAX=1e2;
        addEye2HessianMIN=1e-20;

        maxDirectionError=1e-7;

        if nargin<5
            addEye2HessianU=1e-9;
            addEye2HessianEq=1e-9;
        else
            addEye2HessianU=addEye2Hessian(1);
            if length(addEye2Hessian)<2
                addEye2HessianEq=addEye2Hessian(end);
            else
                addEye2HessianEq=addEye2Hessian(2);
            end
        end
        setAddEye2HessianU__(obj,addEye2HessianU);
        setAddEye2HessianEq__(obj,addEye2HessianEq);
        updateAddEye2HessianU=false;
        updateAddEye2HessianEq=false;
    else
        addEye2HessianU=nan;
        addEye2HessianEq=nan;
    end

    if ~obj.useLDL
        obj.useInertia=false; % cannot use inertia without LDL factorization
    end

    mpDesired=obj.nU;
    if obj.smallerNewtonMatrix
        mnDesired=obj.nG;
    else
        mnDesired=obj.nF+obj.nG;
    end

    %%%%%%%%%%%%%%%%%%
    %% Start timing %%
    %%%%%%%%%%%%%%%%%%

    dt0=clock();
    dt1=clock();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize primal variables & scaling %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize dual variables %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        initDualEqX__(obj);
        %initDualEq__(obj);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Print & initilize headers %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    printf2('%s.m (coupledAlphas=%d,skipAffine=%d,delta=%g,addEye2Hessian=%d,adjustAddEye2Hessian=%d,useInertia=%d,muFactorAggresive=%g,muFactorConservative=%g,useLDL=%d,umfpack=%d):\n   %d primal variable, %d equality constraints, %d inequality constraints\n',...
            FUNCTION__,obj.coupledAlphas,obj.skipAffine,obj.delta,...
            obj.setAddEye2Hessian,obj.adjustAddEye2Hessian,obj.useInertia,...
            obj.muFactorAggressive,obj.muFactorConservative,...
            obj.useLDL,obj.useUmfpack,...
            obj.nU,obj.nG,obj.nF);
    if obj.verboseLevel>=3
        headers='Iter      cost   |grad|   |eq|    ineq.    dual    gap   l(mu) ';
        if obj.setAddEye2Hessian
            headers=[headers,'l(Hp) l(H=) '];
        end
        if obj.adjustAddEye2Hessian
            if obj.useInertia
                headers=[headers,'eig+ eig-  d.err. '];
            else
                headers=[headers,'curvature  d.err. '];
            end
        end
        headers=sprintf('%salphaA  sigma   alphaP  alphaDI alphaDE      time\n',headers);
        if obj.nF>0
            headers=sprintf('%s%4d:<-mx des.->%8.1e%8.1e                %8.1e%6.1f',...
                            headers,maxIter,obj.gradTolerance,obj.equalTolerance,desiredDualityGap,log10(muMin));
        else
            headers=sprintf('%s%4d:<-mx tol.->%8.1e%8.1e                              ',...
                            headers,maxIter,obj.gradTolerance,obj.equalTolerance);
        end
        if obj.setAddEye2Hessian && obj.adjustAddEye2Hessian
            headers=sprintf('%s%6.1f',headers,log10(obj.addEye2HessianUtolerance));
            if obj.useInertia
                headers=sprintf('%s      %5d%5d%8.1e',headers,mpDesired,mnDesired,maxDirectionError);
            else
                headers=sprintf('%s                %8.1e',headers,maxDirectionError);
            end
        end
        fprintf('%s\n',headers);
    end

    %%%%%%%%%%%%%%%
    %% Main loop %%
    %%%%%%%%%%%%%%%

    if obj.debugConvergence
        lastJ=inf;
    end

    iter=0;
    while true
        iter=iter+1;
        if obj.verboseLevel>=3
            if mod(iter,50)==0
                fprintf('%s\n',headers);
            end
            fprintf('%4d:',iter);
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
            J=getJ__(obj);
        end
        if obj.verboseLevel>=3
            fprintf('%11.3e',full(J));
        end

        norminf_grad=getNorminf_Grad__(obj);
        printf3('%8.1e',full(norminf_grad));

        if isnan(norminf_grad)
            printf2('  -> failed to invert hessian\n');
            status = 4;
            break;
        end

        if obj.nG>0
            norminf_eq=getNorminf_G__(obj);
            printf3('%8.1e',full(norminf_eq));
        else
            printf3('  -eq-  ');
        end

        if obj.nF>0
            [gap,ineq,dual]=getGapMinFMinLambda__(obj);
            printf3('%8.1e%8.1e%8.1e',full(ineq),full(dual),full(gap));
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
            printf3(' -ineq-  -dual-   -gap- ');
        end

        if norminf_grad<=obj.gradTolerance && ...
                (obj.nF==0 || gap<=desiredDualityGap) && ...
                (obj.nG==0 || norminf_eq<=obj.equalTolerance) && ...
                (~obj.setAddEye2Hessian || ~obj.adjustAddEye2Hessian || ...
                 addEye2HessianU<=obj.addEye2HessianUtolerance)
            printf2('  -> clean exit\n');
            status = 0;
            break;
        end

        if obj.nF>0
            printf3('%6.1f',log10(mu));
        else
            printf3(' -mu- ');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Adjust addEye2Hessian %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        if obj.setAddEye2Hessian && obj.adjustAddEye2Hessian

            % updates delayed from the last iteration
            if updateAddEye2HessianU
                setAddEye2HessianU__(obj,addEye2HessianU);
                updateAddEye2HessianU=false;
            end
            if updateAddEye2HessianEq
                setAddEye2HessianEq__(obj,addEye2HessianEq);
                updateAddEye2HessianEq=false;
            end

            curvature=getCurvature__(obj); % skip inertia test if curvature looks good (see Zavala and Chiang, 2014)
            if curvature<=0 && obj.useInertia
                [mp,mn]=getHessInertia__(obj);
            end
            derr=getDirectionError__(obj);
            if curvature>0 || (obj.useInertia && mp==mpDesired && mn==mnDesired )
                if ~obj.useInertia
                    printf3('%6.1f%6.1f %+8.1e %8.1e',...
                            log10(addEye2HessianU),log10(addEye2HessianEq),full(curvature),full(derr));
                elseif curvature>0
                    printf3('%6.1f%6.1f c=%+6.0e %8.1e',...
                            log10(addEye2HessianU),log10(addEye2HessianEq),full(curvature),full(derr));
                else
                    printf3('%6.1f%6.1f%5.0f%5.0f%8.1e',...
                            log10(addEye2HessianU),log10(addEye2HessianEq),full(mp),full(mn),full(derr));
                end
                if addEye2HessianU>addEye2HessianMIN && derr<maxDirectionError
                %if addEye2HessianU>addEye2HessianMIN && norminf_grad<=10*obj.gradTolerance
                    addEye2HessianU=max(.75*addEye2HessianU,addEye2HessianMIN);
                    updateAddEye2HessianU=true; % update at next iteration
                end
                if addEye2HessianEq>addEye2HessianMIN && derr<maxDirectionError
                    addEye2HessianEq=max(.75*addEye2HessianEq,addEye2HessianMIN);
                    updateAddEye2HessianEq=true; % update at next iteration
                end
                if addEye2HessianU<addEye2HessianUMAX && derr>maxDirectionError
                    addEye2HessianU=min(10*addEye2HessianU,addEye2HessianUMAX);
                    updateAddEye2HessianU=true; % update at next iteration
                end
                if addEye2HessianEq<addEye2HessianEqMAX && derr>maxDirectionError
                    addEye2HessianEq=min(10*addEye2HessianEq,addEye2HessianEqMAX);
                    updateAddEye2HessianEq=true; % update at next iteration
                end
            else
                for ii=1:1%20
                    change=false;
                    if curvature<0
                        if ~obj.useInertia
                            % inertia not available, increase both
                            if obj.verboseLevel>=4
                                fprintf('%6.1f%6.1f %+8.1e %8.1e\n                                                              ',...
                                        log10(addEye2HessianU),log10(addEye2HessianEq),full(curvature),full(derr));
                            end
                            if addEye2HessianU<addEye2HessianUMAX
                                addEye2HessianU= min(10*max(addEye2HessianU,addEye2HessianMIN),addEye2HessianUMAX);
                                setAddEye2HessianU__(obj,addEye2HessianU);
                                change=true;
                            end
                            if addEye2HessianEq<addEye2HessianEqMAX
                                addEye2HessianEq= min(10*max(addEye2HessianEq,addEye2HessianMIN),addEye2HessianEqMAX);
                                setAddEye2HessianEq__(obj,addEye2HessianEq);
                                change=true;
                            end
                        elseif mp<mpDesired
                            % not enough positive eigenvales
                            if obj.verboseLevel>=4
                                fprintf('%6.1f*%5.1f%5.0f*%4.0f%8.1e\n                                                              ',...
                                        log10(addEye2HessianU),log10(addEye2HessianEq),full(mp),full(mn),full(derr));
                            end
                            if addEye2HessianU<addEye2HessianUMAX
                                addEye2HessianU= min(10*max(addEye2HessianU,addEye2HessianMIN),addEye2HessianUMAX);
                                setAddEye2HessianU__(obj,addEye2HessianU);
                                change=true;
                            end
                            if addEye2HessianEq<addEye2HessianEqMAX
                                addEye2HessianEq= min(2*max(addEye2HessianEq,addEye2HessianMIN),addEye2HessianEqMAX);
                                setAddEye2HessianEq__(obj,addEye2HessianEq);
                                change=true;
                            end
                        elseif mn<mnDesired
                            % not enough negative eigenvales
                            if obj.verboseLevel>=4
                                fprintf('%6.1f%6.1f*%4.0f%5.0f*%7.1e\n                                                              ',...
                                        log10(addEye2HessianU),log10(addEye2HessianEq),full(mp),full(mn),full(derr));
                            end
                            if addEye2HessianU<addEye2HessianUMAX
                                addEye2HessianU= min(2*max(addEye2HessianU,addEye2HessianMIN),addEye2HessianUMAX);
                                setAddEye2HessianU__(obj,addEye2HessianU);
                                change=true;
                            end
                            if addEye2HessianEq<addEye2HessianEqMAX
                                addEye2HessianEq= min(10*max(addEye2HessianEq,addEye2HessianMIN),addEye2HessianEqMAX);
                                setAddEye2HessianEq__(obj,addEye2HessianEq);
                                change=true;
                            end
                        end
                    end
                    if ~change
                        break;
                    end
                    curvature=getCurvature__(obj); % skip inertia test if curvature looks good (see Zavala and Chiang, 2014)
                    if curvature<0 && obj.useInertia
                        [mp,mn]=getHessInertia__(obj);
                    end
                    derr=getDirectionError__(obj);
                end
                if ~obj.useInertia
                    printf3('%6.1f%6.1f %+8.1e %8.1e',...
                            log10(addEye2HessianU),log10(addEye2HessianEq),full(curvature),full(derr));
                elseif curvature>0
                    printf3('%6.1f%6.1f c=%+6.0e %8.1e',...
                            log10(addEye2HessianU),log10(addEye2HessianEq),full(curvature),full(derr));
                else
                    printf3('%6.1f%6.1f%5.0f%5.0f%8.1e',...
                            log10(addEye2HessianU),log10(addEye2HessianEq),full(mp),full(mn),full(derr));
                end
            end
        elseif obj.setAddEye2Hessian
            printf3('%6.1f%6.1f',log10(addEye2HessianU),log10(addEye2HessianEq));
        end

        if obj.debugConvergence
            Lf=getLf__(obj);
            fprintf(' Lf = %9.2e ',Lf);

            %% J large
            J=getJ__(obj);
            if J>obj.debugConvergenceThreshold
                printf2('\n%3d: ATTENTION: cost > %10.2e - cost is probably too large\n',...
                            iter,obj.debugConvergenceThreshold);
            end
            %% Hessian with Large entriesarge
            Hess=getHess__(obj);
            tol=1e6;
            [i1,i2]=find(Hess>tol);
            if ~isempty(i1)
                fprintf('\nATTENTION: Hessian has %d very large entries\n',length(i1));
                for j=1:length(i1)
                    fprintf('\tHess(%4d,%4d)=%8.1e\n',i1(j),i2(j),full(Hess(i1(j),i2(j))));
                end
            end

            %% Hessian singular
            tol=1e-7;
            [vv,eg]=eig(full(Hess),'vector');
            k0=(eg<tol & eg>-tol);
            if any(k0)
                kp=eg>tol;
                km=eg<=-tol;
                fprintf('\nATTENTION: Hessian is singular, with eigenvalues: %4d positive ([%8.1e,%8.1e]), %4d zero ([%8.1e,%8.1e]), %4d negative  ([%8.1e,%8.1e]), expected %d,%d,%d: increase addEye2Hessian\n',...
                        sum(kp),min(eg(kp)),max(eg(kp)),...
                        sum(k0),min(eg(k0)),max(eg(k0)),...
                        sum(km),min(eg(km)),max(eg(km)),...
                        obj.nU,0,obj.nF+obj.nG...
                        );
                if obj.verboseLevel>=5
                    fprintf('       Hessian kernel:\n');
                    disp(vv(:,kk));
                end
                fprintf('              ');
            end
            %% Hessian with wrong trace
            dHess=getdHess__(obj);
            kp=dHess>tol;
            km=dHess<=-tol;
            k0=~kp & ~ km;
            if sum(kp)~=obj.nU || sum(km)~=obj.nF+obj.nG
                fprintf('\nATTENTION: ldl_d(Hessian) has %4d positive entries ([%8.1e,%8.1e]), %4d zero ([%8.1e,%8.1e]), %4d negative  ([%8.1e,%8.1e]), expected %d,%d,%d: increase addEye2Hessian\n',...
                        sum(kp),min(dHess(kp)),max(dHess(kp)),...
                        sum(k0),min(dHess(k0)),max(dHess(k0)),...
                        sum(km),min(dHess(km)),max(dHess(km)),...
                        obj.nU,0,obj.nF+obj.nG...
                        );
                fprintf('              ');
            end
            %% Lf not positive definite
            tol=1e-2;
            Lfuu=getLfuu__(obj);
            [vv,egLfuu]=eig(full(Lfuu),'vector');
            km=egLfuu<=-tol;
            if any(km)
                kp=egLfuu>tol;
                k0=~kp & ~km;
                fprintf('\nATTENTION: Lf not strictly convex, Lfuu eigenvalues: %4d positive ([%8.1e,%8.1e]), %4d zero ([%8.1e,%8.1e]), %4d negative  ([%8.1e,%8.1e])\n',...
                        sum(kp),min(egLfuu(kp)),max(egLfuu(kp)),...
                        sum(k0),min(egLfuu(k0)),max(egLfuu(k0)),...
                        sum(km),min(egLfuu(km)),max(egLfuu(km)));
                if obj.verboseLevel>=5
                    fprintf('       Lfuu negative semidefinite subspace:\n');
                    disp(vv(:,km));
                end
                fprintf('              ');
                %status=55;
                %break;
            end
            % H1=Hess(1:obj.nU,1:obj.nU);
            % H2=Hess(obj.nU+1:end,obj.nU+1:end);
            % H12=Hess(1:obj.nU,obj.nU+1:end);
            % format shorte
            % sort(eig(full(H1))'),;
            % sort(eig(full(H2))'),;
            % sort(eig(full(Hess))'),;
            % format short

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

            g_th=max(abs(J-lastJ)*obj.debugConvergenceThreshold,10*obj.gradTolerance);
            if abs(J-lastJ)<1/obj.debugConvergenceThreshold && norminf_grad>g_th
                grad_=getGrad__(obj);
                k=find(abs(grad_)>g_th);
                if ~isempty(k)
                    printf2('\n%3d: ATTENTION: |Delta J| < %10.2e, but |grad| > %10.2e for %4d entries - poorly conditioned hessian\n',...
                            iter,1/obj.debugConvergenceThreshold,g_th,length(k));
                    for ii=k(:)'
                        printf2('\t grad %3d: %9.2e\n',ii,grad_(ii));
                    end
                end
            end
            lastJ=J;
        end

        if obj.nF==0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% No inequality constraints case %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            setAlphaPrimal__(obj,obj.alphaMax);
            if obj.nG>0
                setAlphaDualEq__(obj,obj.alphaMax);
            end
            printf3('  -alpA- -sigm- ');
            printf3('%8.1e  -nan-    -na-    ',obj.alphaMax);

            updatePrimalDual__(obj);
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Inequality constraints case %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if obj.debugConvergence
                [oldF,oldLambda]=getFLambda__(obj);
                oldmu=mu;
            end
            if obj.skipAffine==1
                printf3(' -alpA-  -sigm- ');
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
                printf3('%8.1e',full(alphaPrimal));

                %% Update mu %%
                
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
                    printf3('%8.1e',sigma);
                    mu=full(max(sigma*gap/obj.nF,muMin));
                    setMu__(obj,mu);
                else
                    printf3(' -sigm- ');
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
                newmu=mu;
            end

            [alphaPrimal,alphaDualIneq]=getMaxAlphas_s__(obj);

            if obj.coupledAlphas && alphaDualIneq<alphaPrimal
                alphaPrimal=alphaDualIneq;
            end

            stepback=.99;
            %stepback=max(1-mu,.9);

            alphaPrimal = stepback * alphaPrimal;

            alphaMax = min(alphaPrimal,obj.alphaMax);

            if (alphaMax >= obj.alphaMin)
                % try max
                alphaPrimal=alphaMax/stepback;
                setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_s__(obj);
                if isnan(ineq)
                    printf2('  -> failed to invert hessian\n');
                    status = 4;
                    break;
                end
                if (ineq>0)
                    % recheck just to be safe in case not convex
                    alphaPrimal = stepback * alphaPrimal;
                    setAlphaPrimal__(obj,alphaPrimal);ineq1=getMinF_s__(obj);
                end
                if ineq<=0 || ineq1<=ineq/10
                    % try min
                    alphaPrimal=obj.alphaMin/stepback;
                    setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_s__(obj);
                    if (ineq>0)
                        % try between min and max
                        alphaPrimal=alphaMax*.95;
                        while alphaPrimal >= obj.alphaMin
                            setAlphaPrimal__(obj,alphaPrimal);ineq=getMinF_s__(obj);
                            if (ineq>0)
                                % backtrace just a little
                                alphaPrimal = stepback * alphaPrimal;
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
                alphaDualIneq = stepback * alphaDualIneq;
                if alphaDualIneq>obj.alphaMax
                    alphaDualIneq = obj.alphaMax;
                end
                %alphaDualEq = obj.alphaMax;
                alphaDualEq = alphaDualIneq;
            end

            if obj.nG>0
                setAlphaDualEq__(obj,alphaDualEq);
            end
            setAlphaDualIneq__(obj,alphaDualIneq);
            updatePrimalDual__(obj);

            if obj.nG>0
                printf3('%8.1e%8.1e%8.1e',full(alphaPrimal),full(alphaDualIneq),full(alphaDualEq));
            else
                printf3('%8.1e%8.1e  -eq-  ',full(alphaPrimal),full(alphaDualIneq));
            end

            %% Update mu %%
                
            if obj.skipAffine==1
                % More aggressive if
                % 1) 'long' newton steps in the affine direction
                % 2) small gradient
                % 3) equality constraints fairly well satisfied
                % (2+3 mean close to the central path)
                %th_grad=norminf_grad<=max(1e-1,1e2*obj.gradTolerance);
                %th_eq=(obj.nG==0) || (norminf_eq<=max(1e-3,1e2*obj.equalTolerance));
                th_grad=            norminf_grad<=max(1e-4,1e0*obj.gradTolerance);
                th_eq=(obj.nG==0) || (norminf_eq<=max(1e-5,1e0*obj.equalTolerance));
                if alphaPrimal>obj.alphaMax/2 && th_grad && th_eq
                    %mu = max(muMin,obj.muFactorAggressive*mu);
                    mu=max(muMin,min(obj.muFactorAggressive*mu,mu^1.5));
                    setMu__(obj,mu);
                    printf3(' * ');
                else
                    if alphaPrimal<.1
                        mu=min(mu0,1.1*mu);
                        setMu__(obj,mu);
                        initDualIneq__(obj);
                        printf3('^');
                    elseif alphaPrimal>.99 && th_eq
                        mu=max(mu*obj.muFactorConservative,muMin);
                        setMu__(obj,mu);
                        printf3('v');
                    else
                        printf3('=');
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
            else
                printf3('   ');
            end

            % if no motion, slowly increase mu
            if (alphaPrimal<obj.alphaMin && alphaDualIneq<obj.alphaMin && alphaDualEq<obj.alphaMin)
                mu=max(mu/obj.muFactorConservative,muMin);
                setMu__(obj,mu);
            end

        end  % if obj.nF==0

        if obj.verboseLevel>=3
            dt1=etime(clock(),dt1);
            fprintf('%8.1fms\n',dt1*1e3);
            dt1=clock();
        end

        %%%%%%%%%%%%%%%%%%%%%%%
        %% Debug small alpha %%
        %%%%%%%%%%%%%%%%%%%%%%%

        if obj.debugConvergence
            if obj.nF>0 && alphaPrimal<obj.alphaMax/5 && isfinite(obj.debugConvergenceThreshold)
                kk=find(newF_s<=0 | newLambda_s<=0);
                fprintf('%3d: ATTENTION: alphaPrimal = %8.2e < %8.2e due to %4d entries (lambda * ineq)\n',...
                        iter,alphaPrimal,obj.alphaMax/5,length(kk));
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

        %%%%%%%%%%%%%%%%%%%%
        %% Check Progress %%
        %%%%%%%%%%%%%%%%%%%%

        if obj.debugConvergence
            [du_,dNu_,dLambda_]=getD__(obj);
            tol=1e5;
            if all(abs(du_./u_)<tol)
                fprintf('%3d: ATTENTION: abs(du/u) < %10.2e for all entries\n',...
                        iter,max(abs(du_./u_)));
            end
            if obj.nG>0
                if all(abs(dNu_./nu_)<tol)
                    fprintf('%3d: ATTENTION: abs(dnu/nu) < %10.2e for all entries\n',...
                            iter,max(abs(dNu_./nu_)));
                end
            end
            if obj.nF>0
                if all(abs(dLambda_./l_)<tol)
                    fprintf('%3d: ATTENTION: abs(dlambda/lambda) < %10.2e for all entries\n',...
                            iter,max(abs(dLambda_./l_)));
                end
            end
        end



    end % while true

    time=etime(clock(),dt0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Complete status output when maxIter reached %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        if (obj.setAddEye2Hessian && obj.adjustAddEye2Hessian & addEye2HessianU>obj.addEye2HessianUtolerance)
                status=bitor(status,2048);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Print exit summary %%
    %%%%%%%%%%%%%%%%%%%%%%%%

    if obj.verboseLevel>=2
        J=getJ__(obj);
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
            if bitand(status,2048)
                fprintf("%clarge addEye2Hessian",sep);
                sep=',';
            end
            fprintf(')\n                ');
        else
            fprintf('%4d:status=0x%s, ',iter,dec2hex(status));
        end
        fprintf('cost=%13.5e, ',full(J));
        norminf_grad=getNorminf_Grad__(obj);
        fprintf('|grad|=%10.2e',full(norminf_grad));
        if obj.setAddEye2Hessian
            fprintf(', addEye2Hessian=[%10.2e,%10.2e]',addEye2HessianU,addEye2HessianEq);
        end
        if obj.nG>0
            fprintf(', |eq|=%10.2e',full(norminf_eq));
        end
        if obj.nF>0
            fprintf(', ineq=%10.2e,\n                 dual=%10.2e, gap=%10.2e, last alpha=%10.2e, last mu=%10.2e',full(ineq),full(dual),full(gap),full(alphaPrimal),mu);
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
