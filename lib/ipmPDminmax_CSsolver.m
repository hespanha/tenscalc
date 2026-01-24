function [varargout]=ipmPDminmax_CSsolver(obj,mu0,maxIter,saveIter,addEye2Hessian)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha, Raphael Chinchilla).  All rights reserved.

FUNCTION__='ipmPDminmax_CSsolver';

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
    addEye2HessianDMAX=1e2;
    addEye2HessianEqMAX=1e2;
    addEye2HessianMIN=1e-20;

    maxDirectionError=1e-7;

    if nargin<5
        addEye2HessianU=1e-9;
        addEye2HessianD=1e-9;
        addEye2HessianEq=1e-9;
    else
        addEye2HessianU=addEye2Hessian(1);
        if length(addEye2Hessian)<2
            addEye2HessianD=addEye2Hessian(end);
        else
            addEye2HessianD=addEye2Hessian(2);
        end
        if length(addEye2Hessian)<3
            addEye2HessianEq=addEye2Hessian(end);
        else
            addEye2HessianEq=addEye2Hessian(3);
        end
    end
    setAddEye2HessianU__(obj,addEye2HessianU);
    setAddEye2HessianD__(obj,addEye2HessianD);
    setAddEye2HessianEq__(obj,addEye2HessianEq);
    updateAddEye2HessianU=false;
    updateAddEye2HessianD=false;
    updateAddEye2HessianEq=false;

    mpUdesired=obj.nU+obj.nGd+obj.nFd;
    mnDdesired=obj.nD;
else
    addEye2HessianU=nan;
    addEye2HessianD=nan;
    addEye2HessianEq=nan;
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
end

if obj.nG>0
    initDualEq__(obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print & initilize headers %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printf2('%s.m (coupledAlphas=%d,addEye2Hessian=%d,adjustAddEye2Hessian=%d,muFactorAggresive=%g,muFactorConservative=%g): %d primal variables (%d+%d), %d eq. constr., %d ineq. constr.\n',...
    FUNCTION__,obj.coupledAlphas,...
    obj.setAddEye2Hessian,obj.adjustAddEye2Hessian,...
    obj.muFactorAggressive,obj.muFactorConservative,...
    obj.nZ,obj.nU,obj.nD,obj.nG,obj.nF);
if obj.verboseLevel>=3
    headers='Iter      cost   |grad|   |eq|    ineq.    dual    gap   l(mu) ';
    if obj.setAddEye2Hessian
        headers=[headers,'l(Hu) l(Hd) l(H=) '];
    end
    if obj.adjustAddEye2Hessian
        headers=[headers,'egU+ egD-  d.err. '];
    end
    headers=sprintf('%salphaP  alphaDI alphaDE      time\n',headers);
    if obj.nF>0
        headers=sprintf('%s%4d:<-mx des.->%8.1e%8.1e                %8.1e%6.1f',...
            headers,maxIter,obj.gradTolerance,obj.equalTolerance,desiredDualityGap,log10(muMin));
    else
        headers=sprintf('%s%4d:<-mx tol.->%8.1e%8.1e                              ',...
            headers,maxIter,obj.gradTolerance,obj.equalTolerance);
    end
    if obj.setAddEye2Hessian && obj.adjustAddEye2Hessian
        headers=sprintf('%s%6.1f%6.1f',headers,...
            log10(obj.addEye2HessianUtolerance),log10(obj.addEye2HessianDtolerance));
        headers=sprintf('%s      %5d%5d%8.1e',headers,mpUdesired,mnDdesired,maxDirectionError);
    end
    fprintf('%s\n',headers);
end

%%%%%%%%%%%%%%%
%% Main loop %%
%%%%%%%%%%%%%%%

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

    if obj.verboseLevel>=3
        f=getf__(obj);
    end
    if obj.verboseLevel>=3
        fprintf('%11.3e',full(f));
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
            (addEye2HessianU<=obj.addEye2HessianUtolerance &&  ...
            addEye2HessianD<=obj.addEye2HessianDtolerance) )
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
        if updateAddEye2HessianD
            setAddEye2HessianD__(obj,addEye2HessianD);
            updateAddEye2HessianD=false;
        end
        if updateAddEye2HessianEq
            setAddEye2HessianEq__(obj,addEye2HessianEq);
            updateAddEye2HessianEq=false;
        end

        for ii=1:30
            derr=getDirectionError__(obj);
            [mpD,mnD]=getHessDinertia__(obj);
            [mpU,mnU]=getHessUinertia__(obj);
            % all looks good?
            if mpU==mpUdesired && mnD==mnDdesired && derr<=maxDirectionError
                if addEye2HessianEq>addEye2HessianMIN
                    addEye2HessianEq=max(.75*addEye2HessianEq,addEye2HessianMIN);
                    updateAddEye2HessianEq=true; % update only at next iteration
                end
                if addEye2HessianU>addEye2HessianMIN
                    addEye2HessianU=max(.75*addEye2HessianU,addEye2HessianMIN);
                    updateAddEye2HessianU=true; % update only at next iteration
                end
                if addEye2HessianD>addEye2HessianMIN
                    addEye2HessianD=max(.75*addEye2HessianD,addEye2HessianMIN);
                    updateAddEye2HessianD=true; % update only at next iteration
                end
                break; % for ii
            else
                change=false;
                % increase updateAddEye2HessianU?
                if mpU<mpUdesired && updateAddEye2HessianU<addEye2HessianUMAX
                    addEye2HessianU= min(10*addEye2HessianU,addEye2HessianUMAX);
                    setAddEye2HessianU__(obj,addEye2HessianU);
                    change=true;
                end
                % increase updateAddEye2HessianD?
                if mnD<mnDdesired && updateAddEye2HessianD<addEye2HessianDMAX
                    addEye2HessianD= min(10*addEye2HessianD,addEye2HessianDMAX);
                    setAddEye2HessianD__(obj,addEye2HessianD);
                    change=true;
                end
                % increase addEye2HessianEq
                if derr>maxDirectionError && addEye2HessianEq<addEye2HessianEqMAX
                    addEye2HessianEq= min(10*addEye2HessianEq,addEye2HessianEqMAX);
                    setAddEye2HessianEq__(obj,addEye2HessianEq);
                    change=true;
                end
                if obj.verboseLevel>=4
                    fprintf('%6.1f%6.1f%6.1f%5.0f%5.0f%8.1e\n                                                              ',...
                        log10(addEye2HessianU),log10(addEye2HessianD),log10(addEye2HessianEq),full(mpU),full(mnD),full(derr));
                end
                if ~change
                    break;
                end
            end
        end
        printf3('%6.1f%6.1f%6.1f%5.0f%5.0f%8.1e',...
            log10(addEye2HessianU),log10(addEye2HessianD),log10(addEye2HessianEq),full(mpU),full(mnD),full(derr));
    else
        if obj.setAddEye2Hessian
            printf3("%6.1f%6.1f%6.1f",...
                log10(addEye2HessianU),log10(addEye2HessianD),log10(addEye2HessianEq));
        end
        if obj.adjustAddEye2Hessian
            derr=getDirectionError__(obj);
            printf3("%8.1e",derr);
        end
    end

    if obj.nF==0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% No inequality constraints case %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setAlphaPrimal__(obj,obj.alphaMax);
        if obj.nG>0
            setAlphaDualEq__(obj,obj.alphaMax);
        end
        printf3('%8.1e  -nan-    -na-    ',obj.alphaMax);

        updatePrimalDual__(obj);
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Inequality constraints case %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [alphaPrimal,alphaDualIneq]=getMaxAlphas_s__(obj);

        if obj.coupledAlphas && alphaDualIneq<alphaPrimal
            alphaPrimal=alphaDualIneq;
        end

        stepback=.99;
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
            %mu=max(muMin,mu*obj.muFactorAggressive);
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
    if (obj.setAddEye2Hessian && obj.adjustAddEye2Hessian & ...
            ( addEye2HessianU>obj.addEye2HessianUtolerance || ...
            addEye2HessianD>obj.addEye2HessianDtolerance ) )
        status=bitor(status,2048);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Print exit summary %%
%%%%%%%%%%%%%%%%%%%%%%%%

if obj.verboseLevel>=2
    f=getf__(obj);
    if obj.nG>0 && status<8 % when status>=8 this has been already been computed
        norminf_eq=getNorminf_G__(obj);
    end
    if obj.nF>0 && status<8 % when status>=8 this has been already been computed
        [gap,ineq,dual]=getGapMinFMinLambda__(obj);
    end

    if (status)
        fprintf('%4d:status=0x%s ',iter,dec2hex(status));
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

    fprintf('cost=%13.5e ',full(f));
    norminf_grad=getNorminf_Grad__(obj);
    fprintf('|grad|=%10.2e',full(norminf_grad));
    if obj.setAddEye2Hessian
        fprintf(', addEye2Hessian=[%10.2e,%10.2e,%10.2e]',...
            addEye2HessianU,addEye2HessianD,addEye2HessianEq);
    end
    if obj.nG>0
        fprintf(', |eq|=%10.2e',full(norminf_eq));
    end
    if obj.nF>0
        fprintf(', ineq=%10.2e,\n                dual=%10.2e, gap=%10.2e, last alpha=%10.2e, last mu=%10.2e',...
            full(ineq),full(dual),full(gap),full(alphaPrimal),mu);
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