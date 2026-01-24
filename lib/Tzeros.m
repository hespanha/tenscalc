function obj=Tzeros(varargin)
% var = Tzeros([])
%
% var = Tzeros([n1,n2,...,na])
%
% var = Tzeros(n1,n2,...,na)
%
%   Returns a Tcalculus tensor with all entries equal to 0.
%   The integers n1,n2,...,na specify the dimension of each index of the tensor.
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
obj=Tcalculus('zeros',osize,[],[],{},1);
end
