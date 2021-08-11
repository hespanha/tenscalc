function obj=myeye(osize)
% var = myeye([])
% var = myeye([n1,n2,...,na,n1,n2,...,na])
%   Returns an identity tensor.  The integers n1,n2,...,na specify the
%   dimension of each index of the tensor.  Note that an identity
%   tensor is expected to have the first half of the dimensions equal
%   to the second half.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if isempty(osize)
        obj=1;
        return;
    end
    ind1=1:length(osize)/2;
    ind2=ind1(end)+1:length(osize);
    if length(ind1)~=length(ind2)
        osize
        ind1,ind2
        error('eye matrix must have an even number of indices');
    end
    if ~isequal(osize(ind1),osize(ind2))
        osize,ind1,ind2
        error('first and second half of dimensions for eye matrix must correspond to compatible sizes');
    end
    obj = zeros(osize);
    if 1
        k=(1:osize(ind1(1)))';
        for i=2:length(ind1)
            k=[kron(ones(osize(ind1(i)),1),k),kron((1:osize(ind1(i)))',ones(size(k,1),1))];
        end
        kk=cell(size(k,2),1);
        for i=1:size(k,2)
            kk{i}=k(:,i);
        end
        ind=sub2ind(osize,kk{:},kk{:});
        obj(ind)=1;
    else
        % untested but should be faster
        osize1=osize(1:length(osize)/2)
        sub1=memory2subscript(osize1,1:prod(osize1))
        [sub1;sub1]
        ind=sub2ind(osize,sub1,sub1)
        error
    end
end
