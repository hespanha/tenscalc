function obj=Tzeros(varargin)
% var = Tzeros([])
%   or 
% var = Tzeros([n1,n2,...,na])
%
%   Returns a Tcalculus tensor with all entries equal to 0.
%   The integers n1,n2,...,na specify the dimension of each index of the tensor.
%
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
    
    if nargin==1
        osize=varargin{1};
    else
        osize=[varargin{:}];
    end
    obj=Tcalculus('zeros',osize,[],[],{},1);
end
