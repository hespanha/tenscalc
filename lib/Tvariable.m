function obj=Tvariable(name,osize,nowarningsamesize,nowarningever)
% var = Tvariable(name,[n1,n2,...,na],nowarningsamesize,nowarningever)
%
% Returns a Tcalculus tensor symbolic variable.  The integers
% n1,n2,...,na specify the dimension of each index of the tensor.
% When the tensor dimensions are missing, they are guessed from the
% value (removing singleton dimensions at the end).
%
% The optional 3rd argument 'nowarningsamesize', prevents a warning
% that is given when one attempts to create a variable that already
% exists. The warning is only omitted when the new variable has the
% same size as the existing one.
%
% The optional 4th argument 'nowarningever', prevents a warning
% that is given when one attempts to create a variable that already
% exists. The warning is always omitted, regardless of the variable sizes.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    if ~ischar(name)
        error('Tvariable: 1st parameter must be a string (variable name)')
    end

    if nargin==1
        % variable name together with size
        s=regexp(name,'(\w+)(\[[^\]]*\])','tokens');
        if length(s)==1 && length(s{1})==2
            osize=s{1}{2};
            name=s{1}{1};
        end
    end

    if ~isvarname(name)
        error('Tvariable: 1st parameter must be a valid variable name (not ''%s'')',name)
    end

    if ~exist('osize','var')
        osize=[];
    end

    if nargin<3
        nowarningsamesize=false;
    end
    if ~islogical(nowarningsamesize)
        error('Tvariable: (optional) 3rd argument must be boolean');
    end

    if nargin<4
        nowarningever=false;
    end
    if ~islogical(nowarningever)
        error('Tvariable: (optional) 4rd argument must be boolean');
    end

    if ischar(osize)
        osize=evalin('caller',osize);
    end

    obj=Tcalculus('variable',osize,name,[],{},1,nowarningsamesize,nowarningever);

    if nargout==0
        assignin('caller',name,obj);
    end
end
