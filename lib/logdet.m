function y=logdet(A)
% logdet - Natural logarithm of the determinant of a matrix
%
%   logdet(A) returns the natural logarithm of the determinant of the
%   square matrix A.
%
    warning('this is a very bad way to cpmpute log-det\n');
    y=log(det(A));
end
