function linear=subscript2memory(tsize,subscript)
% linear=subscript2memory(tsize,subscript)
%    Returns an array of 1-based linear indexing, corresponding to 1-based subscripts,
%    for a tensor of size 'tsize'.
%    Somewhat similar to Matlab's 'sub2ind' but takes a matrix with the subcripts.
% Inputs:
%    tsize (1xn)  - tensor size
%    subscripts (nxN) - matrix with 1-based subscripts, each column corresonding
%                       to one entry of tsize
% Output:
%    linear (1xN) - vector with 1-based indexes into the tensor
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    linear=ones(1,size(subscript,2),'uint64');
    k=uint64([1 cumprod(tsize(1:end-1))]);
    subscript=uint64(subscript);
    for i = length(tsize):-1:1,
        linear=linear+(subscript(i,:)-1)*k(i);
    end
end