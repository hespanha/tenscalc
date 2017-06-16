function obj=Tvariable(name,osize)
% var = Tconstant(name,[n1,n2,...,na])
%
% Returns a Tcalculus tensor symbolic variable.  The integers
% n1,n2,...,na specify the dimension of each index of the tensor.
% When the tensor dimensions are missing, they are guessed from the
% value (removing singleton dimensions at the end).
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
    
    if ~ischar(name)
        error('Tvariable: 1st parameter must be a string (variable name)')
    end
    if nargin<2
        error('Tvariable: missing 2nd parameter (size)');
    end
    
    if ~isvarname(name)
        error('Tvariable: 1st parameter a valid variable name')
    end

    if ischar(osize)
        osize=evalin('caller',osize);
    end

    obj=Tcalculus('variable',osize,name,[],{},1);
    
    assignin('caller',name,obj);
end

