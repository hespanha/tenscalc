classdef fasttable < matlab.mixin.Copyable % abstract handle class with a copy method
% Class that implements a table-like structure optimized for speed of
%  . appending rows
%  . finding rows
%
% Each variable can be
%  . (protected) categorical array
%  . row vector of constant length
%  . cell string
%
% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.

    properties (GetAccess={},SetAccess={})
        data={};          % cell-vector where data will be stored, one variable per entry
        hash=[];          % hash value for type 1,2,3 for faster search (just the sum)
        nRows=0;          % number of rows in table
        dataSize=0;       % number of rows in data arrays
    end

    properties (GetAccess={},SetAccess='immutable')
        variableNames={}; % string array with variable names
        variableTypes=[]; % array with variable types:
                          % 0 - categorical array
                          % 1 - numerical row vector
                          % 2 - cell string
                          % 3 - numeric matrix
                          %     (to be packed into a string before storing)
                          % 4 - general
                          %     (stored unpacked as a cell array)
        variableCategories={};
        nVariables;       % number of variables;

        args={};          % arguments used to create table

        rowBlocks=3000;   % number of rows by which the data arrays are increased
    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method for load and deserialize (must be static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj1=loadobj(obj)
            obj1=fasttable(obj.args{:});
            obj1.data=obj.data;
            obj1.hash=obj.hash;
            obj1.nRows=obj.nRows;
            obj1.dataSize=obj.dataSize;
            %fprintf('loadobj(fasttable): %s -> %s\n',class(obj),class(obj1))
        end
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Create table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=fasttable(varargin)
        % obj=fasttable('var name1',type1,'variable name 2',type2,...)
        %
        %    Creates an table with the given variable names and
        %    types. The types can be:
        %      {cell-string}  = variable will be a categorical array
        %                     with the given categories
        %                     (data stored as a numeric 1-based index)
        %
        %      'fixed-vector' = variable is a row vector always with the same size
        %                     (data stored as numeric matrix)
        %
        %      'string'     = variable is a string
        %                     (data stored as a string array)
        %
        %      'matrix'     = variable is a matrix (different rows may
        %                     have different sizes)
        %                     (packed/unpacked using getByteStreamFromArray/getArrayFromByteStream)
        %
        %      'general'    = variable can be of any type (stored as cell array)
        %                     (packed/unpacked using getByteStreamFromArray/getArrayFromByteStream)
            obj.args=varargin;
            for i=1:2:length(varargin)
                obj.variableNames{end+1,1}=varargin{i};
                if iscell(varargin{i+1})
                    % categorical
                    obj.variableTypes(end+1,1)=0;
                    obj.variableCategories{end+1,1}=varargin{i+1};
                    obj.data{end+1,1}=[];
                elseif isequal(varargin{i+1},'fixed-vector');
                    obj.variableTypes(end+1,1)=1;
                    obj.variableCategories{end+1,1}={};
                    obj.data{end+1,1}=[];
                elseif isequal(varargin{i+1},'string');
                    obj.variableTypes(end+1,1)=2;
                    obj.variableCategories{end+1,1}={};
                    obj.data{end+1,1}={};
                elseif isequal(varargin{i+1},'matrix');
                    obj.variableTypes(end+1,1)=3;
                    obj.variableCategories{end+1,1}={};
                    obj.data{end+1,1}={};
                elseif isequal(varargin{i+1},'general');
                    obj.variableTypes(end+1,1)=3;
                    obj.variableCategories{end+1,1}={};
                    obj.data{end+1,1}={};
                end
            end
            obj.nVariables=length(obj.variableNames);
            obj.hash(end+1,obj.nVariables)=nan;
        end % fasttable

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method for save and serialize
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj1=saveobj(obj)
            obj1=struct('args',{obj.args},...
                        'data',{obj.data},...
                        'hash',{obj.hash},...
                        'nRows',{obj.nRows},...
                        'dataSize',{obj.dataSize});
            %fprintf('saveobj(fasttable): %s -> %s\n',class(obj),class(obj1))
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Access properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function disp(obj,rows)
            if nargin<2
                rows=1:obj.nRows;
            end
            for row=rows
                fprintf('Row %d\n',row);
                for i=1:obj.nVariables
                    switch(obj.variableTypes(i))
                      case 0
                        fprintf('  %s = %g\n',obj.variableNames{i},obj.data{i}(row));
                      case 1
                        fprintf('  %s = %s\n',...
                                obj.variableNames{i},mat2str(obj.data{i}(row,:)));
                      case 2
                        fprintf('  %s = ''%s''\n',...
                                obj.variableNames{i},obj.data{i}{row});
                      case 3
                        %value=str2num(obj.data{i}{row});
                        value=getArrayFromByteStream(uint8(obj.data{i}{row}));
                        if iscell(value)
                            for j=1:length(value)
                                fprintf('  %s{%d} = %sn',...
                                        obj.variableNames{i},j,mat2str(value{j}));
                            end
                        elseif isnumeric(value)
                            fprintf('  %s = %s\n',...
                                    obj.variableNames{i},mat2str(value));
                        else
                            fprintf('  %s =',...
                                    obj.variableNames{i});
                            disp(value)
                        end
                      otherwise
                        fprintf('  %s = ?? (type = %d)\n',obj.variableNames{i},obj.variableTypes(i));
                    end
                end
            end
        end

        function width=width(obj)
            width=obj.nVariables;
        end

        function height=height(obj)
            height=obj.nRows;
        end

        function size=size(obj,dim)
            error('size(fasttable) not implemented: inefficient');
        end


        function rc=isempty(obj)
        % rc=isempty(obj)
        %
        %   Returns true if the table has no rows.
            rc=(obj.nRows==0);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Read table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function value=getOne(obj,variable,row)
        % value=getOne(obj,variable,rows)
        %
        %   Returns the value of a given variable for a given row

            % if length(row) ~= 1
            %     row
            %     error('getOne: a single row expected');
            % end
            i=find(strcmp(variable,obj.variableNames));
            % does not test i for empty for speed
            switch(obj.variableTypes(i))
              case 3
                %value=str2num(obj.data{i}{row});
                value=getArrayFromByteStream(uint8(obj.data{i}{row}));
              case 1
                value=obj.data{i}(row);
              case 0
                value=obj.variableCategories{i}{obj.data{i}(row)};
              case 2
                value=obj.data{i}{row};
            end
        end

        function value=getMulti(obj,variable,rows)
        % value=getMulti(obj,variable,rows)
        %
        %   Returns the value of a given variable for a set of rows
        %   (a little slower than getOne() for a single row)

            if nargin<3
                rows=1:obj.nRows;
            end
            i=find(strcmp(variable,obj.variableNames));
            switch(obj.variableTypes(i))
              case 3
                if isempty(rows)
                    value={};
                else
                    value=obj.data{i}(rows,:);
                    %value=cellfun(@(x)str2num(x),value,'uniform',false);
                    value=cellfun(@(x)getArrayFromByteStream(uint8(x)),value,'uniform',false);
                end
              case 1
                if isempty(rows)
                    value=[];
                else
                    value=obj.data{i}(rows,:);
                end
              case 0
                if isempty(rows)
                    value=[];
                else
                    value=obj.variableCategories{i}(obj.data{i}(rows));
                    value=value(:);
                end
              case 2
                if isempty(rows)
                    value={};
                else
                    value=obj.data{i}(rows,:);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Write to table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function set(obj,variable,row,value)
        % set(obj,variable,row,value)
        %
        %   Sets value of a given variable,, for a given row

            i=find(strcmp(variable,obj.variableNames));
            switch(obj.variableTypes(i))
              case 3
                %value=mymat2str(value);
                value=char(getByteStreamFromArray(value));
                obj.data{i}{row}=value;
                obj.hash(row,i)=sum(value);
              case 1
                obj.data{i}(row)=value;
                obj.hash(row,i)=sum(value);
              case 0
                obj.data{i}(row)=find(strcmp(value,obj.variableCategories{i}));
                obj.hash(row,i)=NaN;
              case 2
                obj.data{i}{row}=value;
                obj.hash(row,i)=sum(value);
            end
        end

        function row=appendRow(obj,varargin)
        % row=appendRow(obj,value1,value2,...)
        %
        %   Appends a row at the end of the table with the given values
        %   for the variables. Returns the index of the row created.

            row=obj.nRows+1;
            obj.nRows=row;
            % reserve memory
            if row>obj.dataSize
                for i=1:length(varargin)
                    switch(obj.variableTypes(i))
                      case 3
                        obj.data{i}(end+1:end+obj.rowBlocks,1)={''};
                      case 1
                        obj.data{i}(end+1:end+obj.rowBlocks,:)=nan(obj.rowBlocks,size(varargin{i},2));
                      case 0
                        obj.data{i}(end+1:end+obj.rowBlocks,1)=nan;
                      case 2
                        obj.data{i}(end+1:end+obj.rowBlocks,1)={''};
                    end
                end
                obj.hash(end+1:end+obj.rowBlocks,obj.nVariables)=nan;
                obj.dataSize=obj.dataSize+obj.rowBlocks;
            end
            variableTypes=obj.variableTypes;
            % set variable values
            for i=1:length(varargin)
                value=varargin{i};
                switch(variableTypes(i))
                  case 3
                    %value=mymat2str(value);
                    value=char(getByteStreamFromArray(value));
                    obj.data{i}{row}=value;
                    obj.hash(row,i)=sum(value);
                  case 1
                    obj.data{i}(row)=value;
                    obj.hash(row,i)=sum(value);
                  case 0
                    obj.data{i}(row)=find(strcmp(value,obj.variableCategories{i}));
                  case 2
                    obj.data{i}{row}=value;
                    obj.hash(row,i)=sum(value);
                end
            end

        end % appendrow

        function row=appendRowUnique(obj,varargin)
        % row=appendRowUnique(obj,value1,match1,value2,match2,...)
        %
        %   Searches for a row that matches the given variables.
        %   When found, returns the row index, otherwise
        %   creates a row with the given variables
        %
        % Inputs:
        %   value1,value2,... - desired values for the variables
        %
        %   match1,match2,... - when true, equality of the variable is
        %                       needed for a match; when false the
        %                       value of the variable is only used if
        %                       a new row is created.
        %
        % Output:
        %   row - index of the row created/updated

            % search for a match
            % handle 1st element separately to go faster
            if varargin{2}
                value=varargin{1};
                %fprintf('checking variable %d (type=%d)\n',i,obj.variableTypes(1));
                switch(obj.variableTypes(1))
                  case 3
                    %value=mymat2str(value);
                    value=char(getByteStreamFromArray(value));
                    row=find(obj.hash(:,1)==sum(value));
                    row=row(strcmp(value,obj.data{1}(row)));
                  case 1
                    row=find(obj.hash(:,1)==sum(value));
                    if ~isempty(row)
                        row=row(all(obj.data{1}(row)==repmat(value,length(row),1),2));
                    end
                  case 0
                    c=find(strcmp(value,obj.variableCategories{1}));
                    if isempty(c)
                        disp(obj.variableCategories{1})
                        error('fasttable: unknown category ''%s''\n',value);
                    end
                    row=find(obj.data{1}==c);
                  case 2
                    row=find(obj.hash(:,1)==sum(value));
                    row=row(strcmp(value,obj.data{1}(row)));
                end
                if isempty(row)
                    row=appendRow(obj,varargin{1:2:end});
                    return
                end
            else
                row=1:obj.nRows;
            end
            for i=2:length(varargin)/2
                if varargin{2*i}
                    value=varargin{2*i-1};
                    %fprintf('checking variable %d (type=%d)\n',i,obj.variableTypes(i));
                    switch(obj.variableTypes(i))
                      case 3
                        %value=mymat2str(value);
                        value=char(getByteStreamFromArray(value));
                        row=row(obj.hash(row,i)==sum(value));
                        row=row(strcmp(value,obj.data{i}(row)));
                      case 1
                        row=row(obj.hash(row,i)==sum(value));
                        row=row(all(obj.data{i}(row)==repmat(value,length(row),1),2));
                      case 0
                        c=find(strcmp(value,obj.variableCategories{i}));
                        if isempty(c)
                            disp(obj.variableCategories{i})
                            error('fasttable: unknown category ''%s''\n',value);
                        end
                        row=row(obj.data{i}(row)==c);
                      case 2
                        row=row(obj.hash(row,i)==sum(value));
                        row=row(strcmp(value,obj.data{i}(row)));
                    end
                    if isempty(row)
                        row=appendRow(obj,varargin{1:2:end});
                        return
                    end
                end
            end

            if length(row)>1
                error('appendRowUnique: multiple matches\n');
            end
        end % appendRowUnique

    end % methods

    methods (Static)
        function t=test()
            clear all
            t=fasttable('var1',{'cat1','cat2'},'var2','string','var3','fixed-vector','var4','matrix')
            N=5000;

            profile on
            t0=clock;
            for i=1:N
                appendRow(t,'cat2','qq',[5,55,7],rand(3,4));
            end
            fprintf('appendRowUnique: %f us/append\n',1e6*etime(clock,t0)/N);
            t

            t0=clock;
            for i=1:N
                appendRowUnique(t,'cat2',true,'tt',true,[5,55,7],true,rand(3,4),false);
            end
            fprintf('appendRowUnique: %f us /append\n',1e6*etime(clock,t0)/N);
            profile viewer
            t

        end %test
    end % methods

end % class
