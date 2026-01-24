function sz=TcheckVariable(name);
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

global TCsymbolicExpressions TCsymbolicExpressionsHash

if ~isempty(TCsymbolicExpressions) && ~strcmp(class(TCsymbolicExpressions),'struct')
    error('unexpected global variable TCsymbolicExpressions\n')
end

types={TCsymbolicExpressions(:).type};
for i=find(strcmp(types,'variable'))
    if myisequal(TCsymbolicExpressions(i).parameters,name)
        sz=TCsymbolicExpressions(i).osize;
        return;
    end
end
sz=[];

end
