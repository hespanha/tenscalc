function obj=norm2(obj1,S)
% norm2 - Squared quadratic norm 
%   
%   norm2(x) returns the sum of the square of all entries of
%   the tensor x, which for matrices corresponds to the square
%   of the Frobenius norm of x.
%   
%   norm2(x,S) returns the value of the quadratic form <x,Sx>. This
%   form is only applicable when x is a vector (tensor with 1
%   dimension) and S a square matrix (tensor with 2 dimensions).

if nargin<2
    obj=sum(obj1(:).^2);
else
    mx=size(obj1);
    ms=size(S);
    if length(mx)==2 && length(ms)==2 && mx(2)==1 && ms(1)==mx(1) && ms(2)==mx(1)
        obj=obj1'*S*obj1;
    else
        error('norm2: norm2(x,S) called with x[%s] not a column vector and/or S[%s] not a square matrix of compatible size\n',index2str(size(obj1)),index2str(size(S)));
    end
end

end

