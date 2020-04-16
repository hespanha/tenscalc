function analyzeHess(Hess,Grad,vars,constr,threshold,u,G,F,nu,lambda,mu)
    
    if nargin<4
        threshold=.1;
    end

    fprintf('Problem has %d variables and %d constraints:\n',length(vars),length(constr));
    primalVariables=struct();
    k=1;
    for i=1:length(vars)
        primalVariables(i).name=name(vars{i});
        primalVariables(i).size=size(vars{i});
        primalVariables(i).range=[k:k+prod(primalVariables(i).size)-1];
        k=primalVariables(i).range(end)+1;
        fprintf('   variable %d: %s[%s]\n',i,primalVariables(i).name,index2str(primalVariables(i).size));
    end
    constraints=struct();
    for i=1:length(constr)
        constraints(i).type=type(constr{i});
        constraints(i).size=size(constr{i});
        constraints(i).range=[k:k+prod(constraints(i).size)-1];
        k=constraints(i).range(end)+1;
        fprintf('   constraint %d: %s [%s]\n',i,constraints(i).type,index2str(constraints(i).size));
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
    for i=1:length(primalVariables)
        ind=find(abs(Grad(primalVariables(i).range))>threshold);
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
    dN=Hess\-[Grad;G;F.*lambda-mu];
    for i=1:length(primalVariables)
        ind=find(abs(dN(primalVariables(i).range))>threshold);
        if ~isempty(ind)
            sub=memory2subscript(primalVariables(i).size,ind);
            fprintf('     Newton direction for variable %s[%s] is large:\n',...
                    primalVariables(i).name,index2str(primalVariables(i).size));
            for q=1:size(sub,2)
                fprintf('        %s(%s)\t= %g\n',...
                        primalVariables(i).name,index2str(sub(:,q)),dN(primalVariables(i).range(ind(q))));
            end
        end
    end
    
    for i=1:length(constraints)
        ind=find(abs(dN(constraints(i).range))>threshold);
        if ~isempty(ind)
            sub=memory2subscript(constraints(i).size,ind);
            fprintf('     Newton direction for multiplier of constraint # %d[%s] of type ''%s'' is large:\n',...
                    i,index2str(constraints(i).size),constraints(i).type);
            for q=1:size(sub,2)
                fprintf('        %s(%s)\t= %g\n',...
                        constraints(i).type,index2str(sub(:,q)),dN(constraints(i).range(ind(q))));
            end
        end
    end

    disp([u,dN(1:length(u))])


end
