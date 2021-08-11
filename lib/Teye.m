function obj=Teye(varargin)
% var = Teye([])
%
% var = Teye([n1,n2,...,na,n1,n2,...,na])
%
% var = Teye(n1,n2,...,na,n1,n2,...,na)
%
%   Returns a Tcalculus identity tensor.  The integers n1,n2,...,na
%   specify the dimension of each index of the tensor.  Note that an
%   identity tensor is expected to have the first half of the
%   dimensions equal to the second half.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if nargin==1
        osize=varargin{1};
    else
        osize=[varargin{:}];
    end
    if mod(length(osize),2)~=0
        osize
        error('eye matrix must have an even number of indices');
    end
    if ~isempty(osize)
        ind1=1:length(osize)/2;
        ind2=ind1(end)+1:length(osize);
    else
        ind1=[];
        ind2=[];
    end
    if ~myisequal(osize(ind1),osize(ind2))
        osize,ind1,ind2
        error('first and second half of dimensions for eye matrix must correspond to compatible sizes');
    end
    if isempty(osize)
        obj =Tcalculus('ones',osize,[],[],{},1);
    else
        obj =Tcalculus('eye',osize,[],[],{},1);
    end
end
