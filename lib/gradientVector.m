function [grad,hess]=gradientVector(objs,vars)
% [grad,hess]=gradientVector(obj,vars)
%
% Return of gradient with respect to a list of variables, reshaped so
% that the derivatives with respect to all variables appear as a
% single additional last dimension (regardless of the size of the
% different variables).
%
% Optionally also return the hessian matrix, also reshaped so that the
% derivatives with respect to all variables appear as two additional
% last dimensions (regardless of the size of the different
% variables).
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if ~iscell(vars)
        vars={vars};
    end
    if ~iscell(objs)
        objs={objs};
    end

    for k=1:length(objs)
        objs{k}=toCalculus(objs{k});
    end

    grad=cell(length(objs),length(vars));
    remove=false(length(vars),1);
    for i=1:length(vars)
        if prod(msize(vars{i}))>0
            if ~isequal(type(vars{i}),'variable')
                vars{i}
                vars{i}.type
                error('%dth variable is not of the type ''variable''\n',i)
            end
            for k=1:length(objs)
                grad{k,i}=gradient(objs{k},vars{i});
            end
        else
            remove(i)=true;
            %fprintf('\ngradientVector: %dth variable ''%s'' is empty\n',i,vars{i}.name);
        end
    end
    if nargout>1
        hess=cell(length(objs),length(vars),length(vars));
        for i=1:length(vars)
            if prod(msize(vars{i}))>0
                for j=i:length(vars)
                    if prod(msize(vars{j}))>0
                        for k=1:length(objs)
                            hess{k,i,j}=gradient(grad{i},vars{j});
                        end
                    end
                end
            end
        end
    end

    vars(remove)=[];
    grad(:,remove)=[];
    if nargout>1
        hess(:,remove,:)=[];
        hess(:,:,remove)=[];
    end

    % "Normalize" shapes
    grad2=cell(size(grad,1),1);
    for k=1:size(grad,1)
        for i=1:size(grad,2)
            % gradient as a vector
            grad{k,i}=reshape(grad{k,i},[size(objs{k}),prod(size(vars{i}))]);
        end
        grad2{k}=cat(length(size(objs{k}))+1,grad{k,:});
    end

    if nargout>1
        hess3=cell(size(grad,1),1);
        for k=1:size(grad,1)
            hess2=cell(size(grad,1),size(grad,2));
            for i=1:size(grad,2)
                for j=i:size(grad,2)
                    % hessian as a 2-matrix
                    hess{k,i,j}=reshape(hess{k,i,j},[size(objs{k}),prod(size(vars{i})),prod(size(vars{j}))]);
                    hess{k,j,i}=hess{k,i,j}';
                end
                hess2{k,i}=cat(length(size(objs{k}))+2,hess{k,i,:});
            end
            hess3{k}=cat(length(size(objs{k}))+1,hess2{k,:});
        end
        hess=cat(length(size(objs{k}))+1,hess2{:});
    end
    grad=cat(1,grad2{:});

end
