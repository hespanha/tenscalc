function outputExpression=packExpressions(inputExpressions)
%  outputExpression=packExpressions(inputExpressions)
%
% Takes as inputs a cell array of Tcalculus expressions and creates a
% single Tcalculus vector (1-index) that is obtained by reshaping all
% the input variables to vectors and stacking them all on top of each
% other.
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

    if isempty(inputExpressions)
        %warning('packExpression: packing an empty input expression\n');
        outputExpression=Tconstant([],[0]);
    else
        for k=1:length(inputExpressions)
            inputExpressions{k}=reshape(inputExpressions{k},prod(size(inputExpressions{k})));
        end
        outputExpression=cat(1,inputExpressions{:});
    end
end