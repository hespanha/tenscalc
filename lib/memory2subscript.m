function subscripts=memory2subscript(tsize,linear)
% subscripts=memory2subscript(tsize,linear)
%    Returns an array of 1-based subscripts, corresponding to 1-based linear indexing
%    for a tensor of size 'tsize'.
%    Somewhat similar to Matlab's 'ind2sub' but returns a matrix with the subcripts.
% Inputs:
%    tsize (1xn)  - tensor size
%    linear (1xN) - vector with 1-based indexes into the tensor
% Output:
%    subscripts (nxN) - matrix with 1-based subscripts, each column corresonding
%                       to one entry of tsize
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

subscripts=zeros(length(tsize),length(linear),'uint64');
k=[1 cumprod(tsize(1:end-1))];
for i = length(tsize):-1:1,
    vi = rem(linear-1, k(i)) + 1;
    subscripts(i,:) = 1+(linear - vi)'/k(i) ;
    linear = vi;
end
end