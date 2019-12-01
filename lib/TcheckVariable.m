function sz=TcheckVariable(name);

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
