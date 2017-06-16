function varargout=toCalculus(varargin)
% [aa,bb,...] = Tconstant(a,b,...)
%
% Converts its argument in TC tensors. If the arguments are already
% TC tensors they remain unchnaged; otherwise they are converted to
% TC tensors using Tconstant.
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

    varargout=varargin;
    for i=1:length(varargin)
        if isa(varargin{i},'Tcalculus')
            continue;
        elseif isnumeric(varargin{i})
            varargout{i}=Tconstant(varargin{i});
            updateFile2table(varargout{i},2);
        else
            varargin{i}
            error('toCalculus: cannot convert class ''%s'' to ''Tcalculus''',class(varargin{i}));
        end
    end
end

