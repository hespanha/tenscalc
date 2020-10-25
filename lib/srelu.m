function y=srelu(x)
% srelu - Soft rectified linear unit activation function.
%
%    srelu(X) returns a tensor with the same size as X, with
%    the "soft" rectified linear unit activation function of
%    the entries of X. Same as log(1+exp(x)) or
%    x+Log(1+exp(-x))
    y=log(1+exp(x));
end