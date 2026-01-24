function [outVariable,whereVariables,packCmd,unpackCmd,varargout]=packVariables(inVariables,outVariableName,varargin)
% [outVariable,whereVariables,packCmd,unpackCmd,out1,out2,... ]=packvariables(inVariables,outVariableName,in1,in2,...)
%
% Takes as inputs a cell array of Tcalculus input variables
% 'inVariables' and creates a single Tcalculus vector (1-index)
% 'outVariable' (named 'outVariableName') that is obtained by
% reshaping all the input variables to vectors and stacking them all
% on top of each other.
%
% The input variables in 'inVariables' are then replaced in all the
% Tcalculus expressions
%    in1,in2,...
% to create new Tcalculus expressions
%    out1,out2,...
% whose dependence on the input variables in 'inVariables' has been
% replaced by approproate dependences on the 'outVariable'.
% See 'help Tcalculus/substitute'
%
% The output 'whereVariables' is a cell array with the indices of
% where each variable is stored in 'outVariable'
%
% The outputs 'packCmd' and 'unpackCmd' return strings with MATLAB
% commands that can be used to pack and unpack a cell array of MATLAB
% numerical variables into the vector.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

nowarningsamesize=true;
nowarningever=true;

verboseLevel=0;

n=0;
packCmd=sprintf('%s=[',outVariableName);
unpackCmd='';
whereVariables=cell(length(inVariables),1);
for i=1:length(inVariables)
    if ~isequal(inVariables{i}.type,'variable')
        inVariables{i}
        inVariables{i}.type
        error('Can only pack variables (not %s)',        inVariables{i}.type);
    end
    len=prod(size(inVariables{i}));
    packCmd=sprintf('%sreshape(%s,%d,1);',packCmd,name(inVariables{i}),len);
    unpackCmd=sprintf('%s%s=reshape(%s(%d:%d),%s);',...
        unpackCmd,inVariables{i}.name,outVariableName,n+1,n+len,...
        index2str(msize(inVariables{i})));
    whereVariables{i}=n+1:n+len;
    n=n+len;
end
packCmd=[packCmd,']'];

outVariable=Tvariable(outVariableName,n,nowarningsamesize,nowarningever);

if verboseLevel>0
    global substituteCounter;
    substituteCounter=0;
end

varargout={};
for i=1:length(varargin)
    varargout{i}=substitute(varargin{i},inVariables,outVariable);
end
if verboseLevel>0
    fprintf('  packVariables: %d substitutions\n',substituteCounter);
end
end
