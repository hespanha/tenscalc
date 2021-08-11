function outputExpression=packExpressions(inputExpressions)
%  outputExpression=packExpressions(inputExpressions)
%
% Takes as inputs a cell array of Tcalculus expressions and creates a
% single Tcalculus vector (1-index) that is obtained by reshaping all
% the input variables to vectors and stacking them all on top of each
% other.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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