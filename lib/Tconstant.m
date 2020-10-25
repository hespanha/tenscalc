function obj=Tconstant(value,osize)
% var = Tconstant(value)
%   or
% var = Tconstant(value,[n1,n2,...,na])
%
% Returns a Tcalculus tensor with a constant value.  The integers
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

    if ~isnumeric(value)
        error('Tconstant: can only accept numeric values (not ''%s'')\n',class(value));
    end

    if nargin<2
        osize=size(value);
        if isequal(osize,[1,1])
            osize=[];
        elseif length(osize)==2 && osize(end)==1
            osize=osize(1);
        end
    end

    % get msize at least 2
    msize=osize;
    while length(msize)<2
        msize(end+1)=1;
    end
    % remove trailing singletons
    while length(msize)>2 && msize(end)==1
        msize(end)=[];
    end

    if ~isequal(msize,size(value)) && (prod(msize)>0 || ~isempty(value))
        disp(value)
        error('Tconstant: size [msize=%s, osize=%s] incompatible with value [%s]\n',index2str(msize),index2str(osize),index2str(size(value)));
    end

    if nnz(value)==0
        obj=Tzeros(osize);
        updateFile2table(obj,1);
    elseif nnz(value)==prod(msize) && isequal(value,ones(msize))
        obj=Tones(osize);
        updateFile2table(obj,1);
    elseif length(osize)==2 && osize(1)==osize(2) && isequal(value,speye(osize))
        obj=Teye(osize);
        updateFile2table(obj,1);
    else
        obj=Tcalculus('constant',osize,value,[],{},1);
    end
end
