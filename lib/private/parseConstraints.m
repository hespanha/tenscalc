function [G,F,nus,lambdas,outputExpressions,template]=parseConstraints(code,classname,constraints,outputExpressions,prefix)
% [G,F,nus,lambdas,template,outputExpressions]=...
%    parseConstraints(code,classname,constraints,outputExpressions,prefix)
%
% Parses an array of constraints into different types in order to
% 1) separates contraints into different type (equality and inequality)
% 2) create the corresponding dual variables
% 3) adds appropriate sets for the dual variables to a csparse
%    object
% 4) generates createGateway templates to generate cmex functions
%    to set values for the dual variable
% 5) replaces appearences of the dual variables in
%    outputExpressions to the given names.
%
%
% Input:
%
% code              - csparse object
% classname         - classname to be used in the names of the functions
% constraints       - cell array of constrainsts as TC symbolic variables
% outputExpressions - original output expressions
% prefix            - string with prefix to be used for the name of
%                     the dual variables
%
% Outputs:
%
% Gcells - cell array of TC symbolic variables corresponding to
%          equality constraints (value should be made equal to
%          zero)
% nus    - cell array of dual variables corresponding to the
%          equality constraints in Gcells
% Fcells - cell array of TC symbolic variables corresponding to
%          inequality constraints (value should be made positive)
% lambdas - cell array of dual variables corresponding to the
%          inequality constraints in Fcells
% outputExpressions - output expressions with the dual variables
%          replaced by the ones that appear in nus & lambdas
% template - struct array with templates to be passed to
%          createGateway() to generate cmex functions to
%          initialize all dual variables (nus and lambdas)
%          (this output argument is optional)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    createSets4Duals=false;

    nowarningsamesize=true;
    nowarningever=true;

    verboseLevel=0;

    if verboseLevel>0
        global substituteCounter;
        substituteCounter=0;
    end


    if nargin<5
        prefix='';
    end

    if nargout>=6
        % template for createGateway
        template=cmextoolsTemplate();
    end

    Gcells={};
    Fcells={};
    nus={};
    lambdas={};
    for k=1:length(constraints)
        if ~isa(constraints{k},'Tcalculus')
            disp(constraints{k})
            error('constraint %d is not a Tcalculus expression\n',k);
        end
        switch (type(constraints{k}))
          case {'iszero'}
            % remove iszero
            op1=Tcalculus(operands(constraints{k}));
            Gcells{end+1}=op1;
            % create appropriate nu
            vname=sprintf('%snu%d_',prefix,length(nus)+1);
            nus{end+1}=Tvariable([vname,'_'],size(op1),nowarningsamesize,nowarningever);
            if nargout>=6 && createSets4Duals
                template(end+1,1).MEXfunction=sprintf('%s_set_%s',classname,name(nus{end}));
                template(end).Cfunction=sprintf('%s_set_%s',classname,name(nus{end}));
                template(end).method=sprintf('setD_%s',name(nus{end}));
                template(end).inputs =struct('type','double',...
                                             'name',name(nus{end}),...
                                             'sizes',size(nus{end}));
            end
            declareSet(code,nus{end},sprintf('%s_set_%s',classname,name(nus{end})));
            if ~isempty(outputExpressions)
                if verboseLevel>1
                    fprintf('substitute: %s -> nus{%d}\n',vname,length(nus));
                end
                outputExpressions=substitute(outputExpressions,...
                                             Tvariable(vname,size(op1),nowarningsamesize,nowarningever),...
                                             nus{end});
            end
          case {'ispositive'}
            % remove ispositive
            op1=Tcalculus(operands(constraints{k}));
            Fcells{end+1}=op1;
            % create appropriate lambda;
            vname=sprintf('%slambda%d_',prefix,length(lambdas)+1);
            lambdas{end+1}=Tvariable([vname,'_'],size(op1),nowarningsamesize,nowarningever);
            if nargout>=6 && createSets4Duals
                template(end+1,1).MEXfunction=sprintf('%s_set_%s',classname,name(lambdas{end}));
                template(end).Cfunction=sprintf('%s_set_%s',classname,name(lambdas{end}));
                template(end).method=sprintf('setD_%s',name(lambdas{end}));
                template(end).inputs =struct('type','double',...
                                             'name',name(lambdas{end}),...
                                             'sizes',size(lambdas{end}));
            end
            declareSet(code,lambdas{end},sprintf('%s_set_%s',classname,name(lambdas{end})));
            if ~isempty(outputExpressions)
                if verboseLevel>1
                    fprintf('substitute: %s -> lambdas{%d}\n',vname,length(lambdas));
                end
                outputExpressions=substitute(outputExpressions,...
                                             Tvariable(vname,size(op1),nowarningsamesize,nowarningever),...
                                             lambdas{end});
            end
          otherwise
            error('constraint of type ''%s'' not implemented\n',type(constraints{k}));
        end
    end

    %% Pack constraints
    G=packExpressions(Gcells);
    F=packExpressions(Fcells);

    if verboseLevel>0
        fprintf('  parseConstraints: %d substitutions\n',substituteCounter);
    end
end
