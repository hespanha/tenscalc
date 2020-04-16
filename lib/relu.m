function y=relu(x)
% relu - Rectified linear unit activation function.
%    
%    relu(X) returns a tensor with the same size as X, with
%    each entry equal the the corresponding entry of X or 0,
%    depending on whether the entry is positive or not. Same
%    as max(X,0).
    y=max(x,0);
end