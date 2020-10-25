function y=traceinv(A)
% traceinv - Trace of the inverse of a matrix
%
%   traceinv(lu(A)) or traceinv(ldl(A)) return the natural
%   logarithm of the determinant of the square matrix A.
    y=trace(inv(A));
end
