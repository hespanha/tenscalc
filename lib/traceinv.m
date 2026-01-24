function y=traceinv(A)
% traceinv - Trace of the inverse of a matrix
%
%   traceinv(lu(A)) or traceinv(ldl(A)) return the natural
%   logarithm of the determinant of the square matrix A.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

y=trace(inv(A));
end
