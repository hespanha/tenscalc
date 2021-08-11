function analyzeHess(Hess,Grad,optimizationVariables,constr,threshold,u,G,F,nu,lambda,mu,dx_s,b_s)
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if nargin<5
        threshold=.1;
    end

    if isstruct(optimizationVariables)
        optimizationVariables=struct2cell(optimizationVariables);
    end

    fprintf('Problem has %d variables and %d constraints:\n',length(optimizationVariables),length(constr));
    primalVariables=struct();
    k=1;
    for i=1:length(optimizationVariables)
        primalVariables(i).name=name(optimizationVariables{i});
        primalVariables(i).size=size(optimizationVariables{i});
        primalVariables(i).range=[k:k+prod(primalVariables(i).size)-1];
        k=primalVariables(i).range(end)+1;
        fprintf('   variable %d:   %-20s\t[%s]\tz[%d..%d]\n',...
            i,primalVariables(i).name,index2str(primalVariables(i).size),...
            primalVariables(i).range(1),primalVariables(i).range(end));
    end
    constraints=struct();
    for i=1:length(constr)
        constraints(i).type=type(constr{i});
        constraints(i).size=size(constr{i});
        if strcmp(constraints(i).type,'iszero');
            % indices for equality constraints
            constraints(i).range=[k:k+prod(constraints(i).size)-1];
            k=constraints(i).range(end)+1;
        end
    end
    for i=1:length(constr)
        if strcmp(constraints(i).type,'ispositive');
            % indices for inequality constraints
            constraints(i).range=[k:k+prod(constraints(i).size)-1];
            k=constraints(i).range(end)+1;
        end
        fprintf('   %s %d:\t\t[%s]\tz[%d..%d]\n',constraints(i).type,i,index2str(constraints(i).size),...
        constraints(i).range(1),constraints(i).range(end));
    end

    fprintf('Analyzing Hessian''s magnitude\n');
    tol=1e6;
    [i1,i2]=find(Hess>tol);
    if ~isempty(i1)
        fprintf('  Hessian has %d very large entries\n',length(i1));
        for j=1:length(i1)
            fprintf('\tHess(%4d,%4d)=%8.1e\n',i1(j),i2(j),full(Hess(i1(j),i2(j))));
        end
    end

    fprintf('Analyzing Hessian''s kernel\n');
    [v,d]=eigs(Hess,4,'smallestabs');
    for j=1:size(v,2)
        fprintf('  Looking at eigenvector %d with eigenvalue = %e\n',j,d(j,j));
        for i=1:length(primalVariables)
            ind=find(abs(v(primalVariables(i).range,j))>threshold);
            if ~isempty(ind)
                sub=memory2subscript(primalVariables(i).size,ind);
                fprintf('     variable %s[%s] has entries along which criteria does not change:\n',...
                        primalVariables(i).name,index2str(primalVariables(i).size));
                for q=1:size(sub,2)
                    fprintf('        %s(%s)\t= %g\n',...
                            primalVariables(i).name,index2str(sub(:,q)),v(primalVariables(i).range(ind(q)),j));
                end
            end
        end

        for i=1:length(constraints)
            ind=find(abs(v(constraints(i).range,j))>threshold);
            if ~isempty(ind)
                sub=memory2subscript(constraints(i).size,ind);
                fprintf('     constraint # %d[%s] of type ''%s'' has entries along which criteria does not change:\n',...
                        i,index2str(constraints(i).size),constraints(i).type);
                for q=1:size(sub,2)
                    fprintf('        %s(%s)\t= %g\n',...
                            constraints(i).type,index2str(sub(:,q)),v(constraints(i).range(ind(q)),j));
                end
            end
        end
    end

    fprintf('Analyzing Gradient\n');
    tol=1e-2;
    for i=1:length(primalVariables)
        ind=find(abs(Grad(primalVariables(i).range))>tol);
        if ~isempty(ind)
            sub=memory2subscript(primalVariables(i).size,ind);
            fprintf('   gradient with respect to variable %s[%s] is large:\n',...
                    primalVariables(i).name,index2str(primalVariables(i).size));
            for q=1:size(sub,2)
                fprintf('      %s(%s)\t= %g\n',...
                        primalVariables(i).name,index2str(sub(:,q)),Grad(primalVariables(i).range(ind(q))));
            end
        end
    end

    fprintf('Analyzing Newton direction\n');
    % for "Large Matrix"
    tol=.25;
    %dN=full(Hess\[-Grad;-G;F-mu./lambda]);
    for i=1:length(primalVariables)
        ind=find(abs(dx_s(primalVariables(i).range)./u(primalVariables(i).range))>tol);
        if ~isempty(ind)
            sub=memory2subscript(primalVariables(i).size,ind);
            fprintf('     Newton direction for variable %s[%s] is large:\n',...
                    primalVariables(i).name,index2str(primalVariables(i).size));
            for q=1:size(sub,2)
                fprintf('        %s(%s)\t= %g + %g\t(%d)\n',...
                        primalVariables(i).name,index2str(sub(:,q)),u(primalVariables(i).range(ind(q))),...
                        dx_s(primalVariables(i).range(ind(q))),primalVariables(i).range(ind(q)));
            end
        end
    end

    x=full([u;nu;lambda]);
    tol=.25;
    for i=1:length(constraints)
        ind=find(abs(dx_s(constraints(i).range)./x(constraints(i).range))>tol);
        if ~isempty(ind)
            sub=memory2subscript(constraints(i).size,ind);
            fprintf('     Newton direction for multiplier of constraint # %d[%s] of type ''%s'' is large:\n',...
                    i,index2str(constraints(i).size),constraints(i).type);
            for q=1:size(sub,2)
                fprintf('        %s(%s)\t= %g + %g\t(%d)\n',...
                        constraints(i).type,index2str(sub(:,q)),x(constraints(i).range(ind(q))),...
                        dx_s(constraints(i).range(ind(q))),constraints(i).range(ind(q)));
            end
        end
    end

    %disp([u,dN(1:length(u))])


end
