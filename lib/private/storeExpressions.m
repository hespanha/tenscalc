%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Storage of intermediate expressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

classdef storeExpressions < handle

    properties
        prefix='';
        functions={};
        expressions={};
        constants={};
    end

    methods

        function obj=storeExpressions(prefix)
            if nargin>0
                obj.prefix=prefix;
            end
        end

        function fname=storeConstant(store,fun)

            for i=1:size(store.constants,1)
                if strcmp(store.constants{i,2},fun)
                    fname=store.constants{i,1};
                    return
                end
            end
            fname=sprintf('%sc%d',store.prefix,size(store.constants,1));
            store.constants(end+1,:)={fname,fun};

        end

        function fname=storeExpression(store,fun)

            for i=1:size(store.expressions,1)
                if strcmp(store.expressions{i,2},fun)
                    fname=store.expressions{i,1};
                    return
                end
            end
            fname=sprintf('%se%d',store.prefix,size(store.expressions,1));
            store.expressions(end+1,:)={fname,fun};

        end

        function [strConstants,strExpressions]=str(store)
            strConstants=[];
            strExpressions=[];
            for i=1:size(store.constants,1)
                strConstants=sprintf('%s  %s=%s\n',strConstants,...
                            store.constants{i,1},store.constants{i,2});
            end
            for i=1:size(store.expressions,1)
                strExpressions=sprintf('%s  %s=%s\n',strExpressions,...
                            store.expressions{i,1},store.expressions{i,2});
            end
        end

        function str=strExpression(store,first)
            str=[];
            if nargin<2
                first=1;
            end
            for i=first:size(store.expressions,1)
                str=sprintf('%s  %s=%s\n',str,...
                            store.expressions{i,1},store.expressions{i,2});
            end
        end

        function disp(store)
            [strConstants,strExpressions]=str(store);
            if ~isempty(strConstants)
                fprintf('constants:\n');
                disp(strConstants)
            end
            if ~isempty(strExpressions)
                fprintf('expressions:\n');
                disp(strExpressions)
            end
        end
    end
end