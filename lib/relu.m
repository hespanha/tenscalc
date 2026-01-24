function y=relu(x)
% relu - Rectified linear unit activation function.
%
%    relu(X) returns a tensor with the same size as X, with
%    each entry equal the the corresponding entry of X or 0,
%    depending on whether the entry is positive or not. Same
%    as max(X,0).
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

y=max(x,0);
end