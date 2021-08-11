function varargout=toCalculus(varargin)
% [aa,bb,...] = Tconstant(a,b,...)
%
% Converts its argument in TC tensors. If the arguments are already
% TC tensors they remain unchnaged; otherwise they are converted to
% TC tensors using Tconstant.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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

