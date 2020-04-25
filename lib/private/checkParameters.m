function parameters=checkParameters(parameters)
    
    if isa(parameters,'Tcalculus')
        parameters={parameters};
    end
    
    if ~iscell(parameters)
        parameters
        error('parameters must be a cell array of Tcalculus variables');
    end
    
    for i=1:length(parameters)
        if ~isequal(class(parameters{i}),'Tcalculus')
            parameters{i}
            error('all parameters must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(parameters{i}));
        end
        if ~isequal(type(parameters{i}),'variable')
            parameters{i}
            error('all parameters must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(parameters{i}));
        end
    end
    
end