function y=logdet(A)
% logdet - Natural logarithm of the determinant of a matrix
%
%   logdet(A) returns the natural logarithm of the determinant of the
%   square matrix A.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

warning('this is a very bad way to cpmpute log-det\n');
y=log(det(A));
end
