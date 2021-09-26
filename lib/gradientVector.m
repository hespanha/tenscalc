function [grad,hess]=gradientVector(objs,vars,dependencies)
% [grad,hess]=gradientVector(obj,vars)
%
% Return of gradient of a list of expressions with respect to a list
% of variables, reshaped so that
%
% 1) all entries of all the expressions in the list appear as a single
%    column vector (regardless of their original sizes)
%
% 2) the derivatives with respect to all variables appear as a single
%    additional last dimension (regardless of the size of the
%    different variables).
%
% Optionally also return the hessian matrix, also reshaped so that the
% derivatives with respect to all variables appear as two additional
% last dimensions (regardless of the size of the different
% variables).
%
% [grad,hess]=gradientVector(obj,vars,dependencies)
%
% The optional cell array 'dependecies' has the same size as objs and each entry
%      dependencies{k} = indices of the variables that affect objs{k}
%
% All variables not listed in dependencies{k} are assumed no to affect
% the value of objs{k} and the correspondiong derivative will be
% returned as zero, without checking. This can greatly increazse
% speed, but will lead to errors if variables are missing from objs{k}
%
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

    if numel(objs)==1 && isempty(size(objs{1}))
        scalar=true;
    else
        scalar=false;
    end

    if nargin>=3 && numel(objs)~=numel(dependencies)
        objs,dependencies,
        error('number of expressions %d, must match number of arrays in dependencies variable %d\n',...
              numel(objs),numel(dependencies));
    end

    gradCell=cell(length(objs),length(vars));
    remove=false(length(vars),1);
    for i=1:length(vars)
        if prod(msize(vars{i}))>0
            if ~isequal(type(vars{i}),'variable')
                vars{i}
                vars{i}.type
                error('%dth variable is not of the type ''variable''\n',i)
            end
            for k=1:length(objs)
                if nargin<3 || ( nargin >= 3 && ismember(i,dependencies{k}))
                    gradCell{k,i}=gradient(objs{k},vars{i});
                else
                    gradCell{k,i}=Tzeros([size(objs{k}),size(vars{i})]);
                end
            end
        else
            remove(i)=true;
            %fprintf('\ngradientVector: %dth variable ''%s'' is empty\n',i,vars{i}.name);
        end
    end
    if nargout>1
        hessCell=cell(length(objs),length(vars),length(vars));
        for i=1:length(vars)
            if prod(msize(vars{i}))>0
                for j=i:length(vars)
                    if prod(msize(vars{j}))>0
                        for k=1:length(objs)
                            if nargin<3 || ( nargin >= 3 && ismember(i,dependencies{k}) && ismember(j,dependencies{k}) )
                                hessCell{k,i,j}=gradient(gradCell{k,i},vars{j});
                            else
                                hessCell{k,i,j}=Tzeros([size(objs{k}),size(vars{i}),size(vars{j})]);
                            end
                        end
                    end
                end
            end
        end
    end

    vars(remove)=[];
    gradCell(:,remove)=[];
    if nargout>1
        hessCell(:,remove,:)=[];
        hessCell(:,:,remove)=[];
    end

    % "Normalize" shapes
    grad2=cell(size(gradCell,1),1);
    for k=1:size(gradCell,1)
        for i=1:size(gradCell,2)
            % gradient as a vector
            if scalar
                newsize=[prod(size(vars{i}))];
            else
                newsize=[numel(objs{k}),prod(size(vars{i}))];
            end
            gradCell{k,i}=reshape(gradCell{k,i},newsize);
        end
        if scalar
            grad2{k}=cat(1,gradCell{k,:});
        else
            grad2{k}=cat(2,gradCell{k,:});
        end
    end
    grad=cat(1,grad2{:});

    if nargout>1
        hess3=cell(size(gradCell,1),1);
        for k=1:size(gradCell,1)
            hess2=cell(size(gradCell,1),size(gradCell,2));
            for i=1:size(gradCell,2)
                for j=i:size(gradCell,2)
                    % hessian as a 2-matrix
                    if scalar
                        newsize=[prod(size(vars{i})),prod(size(vars{j}))];
                    else
                        newsize=[numel(objs{k}),prod(size(vars{i})),prod(size(vars{j}))];
                    end
                    hessCell{k,i,j}=reshape(hessCell{k,i,j},newsize);
                    % get hessCell{k,j,i} by transpose
                    if numel(newsize)==2
                        hessCell{k,j,i}=hessCell{k,i,j}';
                    else
                        hessCell{k,j,i}=tprod(hessCell{k,i,j},[1,3,2]);
                    end
                end
                if scalar
                    hess2{k,i}=cat(2,hessCell{k,i,:});
                else
                    hess2{k,i}=cat(3,hessCell{k,i,:});
                end
            end
            if scalar
                hess3{k}=cat(1,hess2{k,:});
            else
                hess3{k}=cat(2,hess2{k,:});
            end
        end
        hess=cat(1,hess3{:});
    end

end
