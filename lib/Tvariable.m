function obj=Tvariable(name,osize,nowarningsamesize)
% var = Tconstant(name,[n1,n2,...,na],nowarningsamesize))
%
% Returns a Tcalculus tensor symbolic variable.  The integers
% n1,n2,...,na specify the dimension of each index of the tensor.
% When the tensor dimensions are missing, they are guessed from the
% value (removing singleton dimensions at the end).
%   
% The optional 3rd argument 'nowarningsamesize', prevents a warning
% that is given when one attempts to create a variable that already
% exists. The warning is only omitted when the new variable has the
% same size as the existing one.
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
        osize=[];
    end
    
    if nargin<3
        nowarningsamesize=false;
    end
    if ~islogical(nowarningsamesize)
        error('Tvariable: (optional) 3rd argument must be boolean');
    end
    
    if ~isvarname(name)
        error('Tvariable: 1st parameter a valid variable name')
    end

    if ischar(osize)
        osize=evalin('caller',osize);
    end

    obj=Tcalculus('variable',osize,name,[],{},1,nowarningsamesize);
    
    assignin('caller',name,obj);
end

