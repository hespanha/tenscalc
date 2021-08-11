function parameters=checkParameters(parameters)
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if isa(parameters,'Tcalculus')
        parameters={parameters};
    end

    if isa(parameters,'struct')
        parameters=struct2cell(parameters);
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