function obj=Tones(varargin)
% var = Tones([])
%
% var = Tones([n1,n2,...,na])
%
% var = Tones(n1,n2,...,na)
%
%   Returns a Tcalculus tensor with all entries equal to 1.
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
    obj=Tcalculus('ones',osize,[],[],{},1);
end
