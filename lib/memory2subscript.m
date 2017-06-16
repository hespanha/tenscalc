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
% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.
    subscripts=zeros(length(tsize),length(linear),'uint64');
    k=[1 cumprod(tsize(1:end-1))];
    for i = length(tsize):-1:1,
        vi = rem(linear-1, k(i)) + 1;
        subscripts(i,:) = 1+(linear - vi)'/k(i) ;
        linear = vi;
    end
end