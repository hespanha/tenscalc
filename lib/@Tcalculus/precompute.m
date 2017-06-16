function varargout=precompute(varargin)
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

    maxNumel2ExpandConstants=2000;
    maxNnz2ExpandConstants=1000;
                             
    varargout=varargin;
    
    for i=1:length(varargin)
        type=getType(varargin{i});
        if ~isequal(type,'constant') && ~isequal(type,'zeros') && ~isequal(type,'ones') && getIsStatic(varargin{i})
            osize=size(varargin{i});
            if prod(osize)>maxNumel2ExpandConstants
                fprintf('precompute: not pre-computing constant [%s] with two many elements (numel = %d>%d)\n',index2str(osize),prod(osize),maxNumel2ExpandConstants);
                continue;
            end
            x=optimize4matlab(varargin{i});
            ex=eval(str(x));
            if 0%length(size(varargin{i}))>2
                fprintf('precompute: not pre-computing constant [%s] with dimension>=2 (numel = %d, nnz = %d)\n',index2str(size(x)),prod(size(x)),nnz(ex));
                continue;
            end
            if nnz(ex)>maxNnz2ExpandConstants
                fprintf('precompute: not precomputing constant [%s] with too many non-zero elements (numel = %d, nnz = %d>%d)\n',index2str(size(x)),prod(size(x)),nnz(ex),maxNnz2ExpandConstants);
                continue;
            end
            %varargin{i}
            varargout{i}=Tconstant(ex,osize);
            fprintf('precompute: precomputing (osize: %s -> %s, nnz: %d)\n',index2str(osize),index2str(size(varargout{i})),nnz(ex));
        end
    end
end

