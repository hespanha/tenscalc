function y=vec2tensor(obj1,sz,subs,dim)
% vec2tensor - Expands a vector to a sparse tensor
%   
%   vec2tensor(X,sz,subs), given 
%     * Tcalculus n-vector X (tensor with size [n]) a
%     * vector sz with d integers
%     * nxd matrix subs of subscripts
%   returns 
%     * Tcalculus tensor Y with size sz, with the nonzero entries
%       taken from X, with
%          Y(subs(i,:))=X(i) for i=1:n
%
%   vec2tensor(X,sz,subs,dim), given 
%     * Tcalculus tensor X with size(X,dim)=n
%     * vector sz with d integers
%     * nxd matrix subs of subscripts
%   returns 
%     * Tcalculus tensor Y with size similar to that of X, but the dim
%       dimension expanded to the sizes in sz, and the nonzero entries
%       taken from X, with
%          Y(...,subs(i,:),...)=X(...,i,...) for i=1:n
%       where the ... denote indices of the dimensions before and
%       after dim
    
    if ndims(obj1)<1
        error('vec2tensor can only be used for tensors with 1 or more dimensions')
    end
    
    if length(dim)~=1
        disp(dim)
        error('vec2tensorc''s 3rd argument must be a scalar')
    end
    
    if dim>ndims(obj1)
        error('vec2tensor''s 3rd argument subs (%d) must be no larger than number of tensor dimensions [%s]\n',...
              dim,index2str(size(obj1)))
    end
    
    if length(sz)~=size(subs,2)
        error('vec2tensor: length of 2nd argument sz (%d) must match number of columns of 3rd argument subs ([%s])',...
              length(sz),index2str(subs));
    end
    
    osize1=size(obj1);
    if size(subs,1)~=osize1(dim)
        error('vec2tensor: number of rows of 3rd argument subs (size(subs)=[%s]) must match size(X[%s],dim=%d)',...
              index2str(size(subs)),index2str(osize1),dim);
    end
    
    for i=1:length(sz)
        if any(subs(:,i)<1)
            error('subsref: subscript in dimension %d smaller than 1\n',i);
        end
        if any(subs(:,i)>sz(i))
            error('vec2tensor: subscript in dimension %d exceeds tensor dimension (%d)\n',i,sz(i))
        end
    end

    % convert osize1 to tenscal sizes
    if osize1(end)==1
        osize1(end)=[];
    end
    
    osize=[osize1(1:dim-1),sz(:)',osize1(dim+1:end)];
    k=find(obj1);
    v=obj1(k);
    subscripts=memory2subscript(osize1,k);
    subscripts=[subscripts(1:dim-1,:);
                subs(subscripts(dim,:),:)';
                subscripts(dim+1:end,:)];
    if size(subscripts,1)==2
        y=sparse(double(subscripts(1,:)),double(subscripts(2,:)),v,osize(1),osize(2));
    else
        k=subscript2memory(osize,subscripts);
        y=zeros(osize);
        y(k)=v;  
    end
end