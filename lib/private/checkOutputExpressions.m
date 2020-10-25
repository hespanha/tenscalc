function [outputExpressions,outputNames]=checkOutputExpressions(outputExpressions)

    if ~iscell(outputExpressions) && ~isstruct(outputExpressions)
        outputExpressions
        error('outputExpressions must be a cell array of Tcalculus variables');
    end

    if iscell(outputExpressions)
        outputNames=cell(length(outputExpressions),1);
        for i=1:length(outputExpressions)
            outputNames{i}=sprintf('y%d',i);
            if ~ismember(class(outputExpressions{i}),{'Tcalculus','double'})
                outputExpressions{i}
                error('outputExpression{%d} is not a Tcalculus variable',i);
            end
        end
    else
        outputNames=fields(outputExpressions);
        for i=1:length(outputNames)
            if ~ismember(class(outputExpressions.(outputNames{i})),{'Tcalculus','double'})
                outputExpressions.(outputNames{i})
                error('outputExpression.%s is not a Tcalculus variable',outputNames{i});
            end
        end
        outputExpressions=struct2cell(outputExpressions);
    end
    [outputExpressions{:}]=toCalculus(outputExpressions{:});
    outputExpressions=outputExpressions(:);
    outputNames=outputNames(:);
end