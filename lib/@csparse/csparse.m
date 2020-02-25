classdef csparse < handle
% Class used to generate C code to perform operations on multi-dimensional matrices.
% The operations to be performed are declared off-line using TCalculus (symbolic) objects.
% 
% See csparse.pdf for details
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
    
    properties
        
        debug; % When nonzero the C functions print debug information
               % to stderr:
               %    >=1 - functions write their name when called
               %         
               %    >=2 - functions write the content of the buffer
               %    >=3 - a get function is declared for every 
               %          intermediate vectorized operations

        tprod2matlab=true; % attemps to convert any tprod to matlab operations
        
        fastRedundancyCheck=true; % uses fast newInstruction some scalarizations

        
        %% Information about set's
        sets=struct('functionName',{},... % desired name for the C function
                    'destination',{},...  % elementary expression to be set
                    'childrenGroups',{}); % groups of instructions that depend on destination
        
        %% Information about get's
        gets=struct('functionName',{},... % desired name for the C function
                    'source',{},...       % elementary expression to be retrieved
                    'parentGroups',{});   % groups of instructions that source depends on
        
        %% Information about copy's
        copies=struct('functionName',{},...   % desired name for the C function
                      'source',{},...         % source elementary expression(s)
                      'destination',{},...    % destination elementary expression(s)
                      'childrenGroups',{},... % groups of instructions that depend on destination
                      'parentGroups',{});     % groups of instructions that source depends on
        
        %% Information about save's
        saves=struct('functionName',{},...  % desired name for the C function
                     'filename',{},...      % filename where variable will be saved
                     'source',{},...        % elementary expression to be retrieved
                     'parentGroups',{},...  % groups of instructions that source depends on
                     'magic',{});           % magic number to match files

        %% Information about ExternalFunctions
        externalFunctions=struct('fileName',{},...  % filename name where source resides
                          'functionName',{},...  % function name
                          'inputs',{},...        % input parameters
                          'output',{},...        % output parameters
                          'defines',{});         % #defines
        
        %% template for createGateway
        template=struct('MEXfunction',{},...% string
                        'Sfunction',{},...  % string
                        'Cfunction',{},...  % string
                        'method',{},...     % string
                        'help',{},...       % string
                        'inputs',struct(...  
                            'type',{},...   % string
                            'name',{},...   % cell-array of strings (one per dimension)
                            'sizes',{}),... % cell-array of strings (one per dimension)
                        'outputs',struct(...% string
                            'type',{},...   % string
                            'name',{},...   % cell-array of strings (one per dimension)
                            'sizes',{}),... % cell-array of strings (one per dimension)
                        'preprocess',{},... % strings (starting with parameters in parenthesis)'
                        'includes',{});     % cell-array of strings (one per file)
        
        %% Information about all the vectorized operations needed:
        %    each vectorized operation is essentially one TC operation.
        
        vectorizedOperations=[];
        nAddedVectorizedOperations=0; % counter of added vectorizedOperations (to keep track of reuse)
        TCindex2CSvectorized=[]; % index in vectorizedOperations of
                                 % entries of TCsymbolicExpressions
        
        %% Information about all the scalar operations needed
        %    each vectorized operation should map to just a few
        %    assembly commands (so that one can reuse instructions
        %    as much as possible).
        % The table is storied in a global variable of the instructionsTable.c
        
        Itypes; % instruction types;
        
        atomicVariables=struct(...
            'Ap',{},...   % column compressed form of indices
            'Ai',{}...
            );
        
        %% Information about dependencies between instructions:
        %    dependencyGraph(i,j) = true if instructions(j) affects instructions(i) 
           dependencyGraph=sparse([],[],[],0,0); 
        % array with the group of each instruction 
        % (instructions in the same group can be computed in blocks)
           instructionsGroup;    
        % dependencyGroup(g,:) = instructions in g dependency group
           dependencyGroups=[]; 
           dependencyGroupColName={};
        % row vector with the memory location for the result of each instruction
           memoryLocations=[]; 
        % list of instructions that appear in sets or copies (as destinations)
           inputInstructions=zeros(1,0);
        % list of instructions that appear in gets or copies (as sources)
           outputInstructions=zeros(1,0);
        % true if instruction appears in gets, sets, copies, or communicate between groups.
           nonLocalInstructions=zeros(1,0);
        % groupLocalInstructions{g} = list of instructions that are only
        %                             used within group g
           groupLocalInstructions={};
        % groupInputInstructions{g} = list of instructions of group g that
        %                             appear in sets or copies (as destinations)
           groupInputInstructions={};
        % groupOutputInstructions{g} = list of instructions of group g that
        %                             appear in gets or copies (as sources)
           groupOutputInstructions={};
        
        %% Compilation information
        
        % C-type for the inermediate computations
        scratchbookType='double'; 
        
        %% Statistics about code generation
        statistics=struct(...
            'nGroups',nan,...          % # of dependency groups  
            'nSets',nan,...            % # of set functions
            'nGets',nan,...            % # of get functions
            'nCopies',nan,...          % # of copy functions
            'nInstructions',nan,...    % final # of instructions
            'nAddedInstructions',0,... % number of (nonunique) instructions added
            'sizeScratchbook',nan,...
            'lu',{{}},...              % information about lu-factorization
            'ldl',{{}},...             % information about ldl-factorization
            'chol',{{}}...             % information about chol-factorization
            );
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method for load and deserialize (must be static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj1=loadobj(obj)
            obj1=csparse(obj.scratchbookType,obj.debug,...
                         obj.tprod2matlab,obj.fastRedundancyCheck);
            fs=fields(obj);
            for i=1:length(fs)
                if strcmp(fs{i},'vectorizedOperations')
                    obj1.vectorizedOperations=fasttable.loadobj(obj.vectorizedOperations);
                else
                    obj1.(fs{i})=obj.(fs{i});
                end
            end
            fprintf('loadobj(csparse): %s -> %s (ATTENTION: Instructions table was lost)\n',class(obj),class(obj1));
        end
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Object creation & display method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=csparse(scratchbookType,debug,...
                             tprod2matlab,fastRedundancyCheck)
            % obj=csparse(scratchbookType,debug,...
            %             tprod2matlab,fastRedundancyCheck)
            %
            % Creates a csparse object. All input parameters are optional.
            %
            % scratchbookType - type of variable used for the
            %                   scratchbook generally 'double', but
            %                   could be 'float'. ATTENTION: float often
            %                   lead to numerical problems.  Default =
            %                   'double'
            %
            % debug - when nonzero, manipulations of the csparse
            %         object will print debug information.
            %         Default = 0
            %
            % tprod2matlab - when true, an attempt is always made to
            %                convert tprod operation to matlab
            %                operations, instead of implementing
            %                tprod directly. Convertion to matlab
            %                sometimes force convertion of sparse
            %                matrices to full.
            %                Default = false
            %
            % fastRedundancyCheck - when true, very intensive
            %                operations (like the lu factorization)
            %                do not check if it is possible to
            %                reuse computations.
            %                Default = false
            
            if nargin<1
                scratchbookType='double';
            end
            if nargin<2
                debug=0;
            end
            if nargin<3
                tprod2matlab=false;
            end
            if nargin<4
                fastRedundancyCheck=false;
            end
            
                
            obj.scratchbookType=scratchbookType;
            obj.debug=debug;
            obj.tprod2matlab=tprod2matlab;
            obj.fastRedundancyCheck=fastRedundancyCheck;

            global TCsymbolicExpressions
            obj.TCindex2CSvectorized=zeros(length(TCsymbolicExpressions),1);
            
            obj.vectorizedOperations=fasttable(...
                'type',{...      % type of operation
                    'variable','zeros','constant','eye','ones',...
                    'subsref','cat','reshape','repmat','full',...
                    'plus','tprod','norm2','norm1','norminf','all','any','min','max','min2','max2',...
                    'abs','clp','compose',...
                    'lu','lu_sym','ldl','ldl_l','ldl_d','chol',...
                    'inv_ldl','logdet_ldl','traceinv_ldl',...
                    'inv_lu','logdet_lu','traceinv_lu',...
                    'mldivide','mldivide_l1','mldivide_u','mldivide_u1','mldivide_d',...
                    'rdivide','mtimes','ctranspose','times','sum','diag','tprod_matlab'...
                       },...
                'description','string',... % description of the operation result
                'name','string',... % (unique) name of the expression:
                ...                 % when type=='variable' is the
                ...                 % name of the variable otherwise,
                ...                 % name is contructed from type
                'osize','matrix',... % size of the sparse array (row
                ...                  % vector)
                'subscripts','general',... % subscripts of the nonzero
                ...       % elements, as a [n x nnz] matrices, where
                ...       %   n   - dimension of the tensor
                ...       %   nnz - number of non-zero entries
                ...       %   each column denotes the index of a nonzero entry                
                'operands','matrix',...   % array with the indices of the operands
                'parameters','general',... % parameters for the instruction
                ...       % . value for type='constant'
                ...       % . array with +/-1 for type='plus'
                ...       % . cell-array with size of sums, followed by indices for type='tprod' 
                ...       % . {row permutation vector (p),column permutation vector (q),...
                ...       %    typical subscripts filename,typical value filename}
                ...       %   for type 'lu'
                ...       % . {row/column permutation vector (p),...
                ...       %    typical subscripts filename,typical value filename}
                ...       %   for type 'ldl' & 'chol'
                ...       % . S structure for subsref
                ...       % . dimension for type='cat, 'min', 'max', 'all'
                'instructions','general',... % [nnz x 1] array with the
                ...          % instructions that compute the nonzero
                ...          % elements referenced in 'subscripts
                ...          % (indices into the 'instructions'
                ...          % table).;
                'atomic','fixed-vector'... % boolean variable indicating
                ...                        % if the operator is atomic
            );
                
            obj.Itypes=instructionTypes();
            
            if ~libisloaded('instructionsTable')
                [notfound,warnings]=loadlibrary('instructionsTable','instructionsTable.h');
            end
            
            initInstructionsTable();
        end
        
        function disp(obj,verboseLevel)
            if nargin<2
                verboseLevel=0;
            end
            fprintf('%s',str(obj,verboseLevel));
        end

        function obj1=copy(obj)
           % Deep copy
            error('copy incomplete\n');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method for save and serialize 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj1=saveobj(obj)
            obj1=struct();
            fs=fields(obj);
            for i=1:length(fs)
                if strcmp(fs{i},'vectorizedOperations')
                    obj1.vectorizedOperations=saveobj(obj.vectorizedOperations);
                else
                    obj1.(fs{i})=obj.(fs{i});
                end
            end
            %obj1
            %fprintf('saveobj(csparse): %s -> %s\n',class(obj),class(obj1))
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Declare operation methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function declareSet(obj,TCdestination,functionName,helpmsg)
        % declareSet(obj,TCdestination,functionName)
        %
        % Declares a 'set' operation to be compiled.  Each 'set'
        % operation loads data from a (full) array into a compiled
        % sparse variable.
            if isempty(TCdestination)
                TCdestination
                error('empty declareSet')
            end
            before=length(obj.vectorizedOperations);
            
            k=addTCexpression(obj,TCdestination);
            obj.sets(end+1,1)=...
                struct('functionName',{functionName},'destination',{k},'childrenGroups',{[]});
            if ~strcmp(type(TCdestination),'variable')
                fprintf('declareSet: ATTENTION destination is not a Tvariable, which may lead to unpredictable behavior\n');
            end
            if length(obj.vectorizedOperations)>before
                fprintf('declareSet: ATTENTION setting a previously unseen variable, which may lead to unpredictable behavior\n');
            end
            
            % add to template
            obj.template(end+1,1).MEXfunction=sprintf('%s',functionName);
            obj.template(end,1).Sfunction=sprintf('%sS',functionName);
            obj.template(end,1).Cfunction=sprintf('%s',functionName);
            obj.template(end,1).method=sprintf('%s',functionName);
            if nargin<4
                obj.template(end,1).help='';
            else
                obj.template(end,1).help=helpmsg;
            end            
            obj.template(end,1).inputs(1).type=obj.scratchbookType;
            if strcmp(type(TCdestination),'variable')
                obj.template(end,1).inputs(1).name=name(TCdestination);
            else
                obj.template(end,1).inputs(1).name='input1';
            end
            obj.template(end,1).inputs(1).sizes=TCdestination.osize;
            
            %computeScalarInstructions(obj,k);
        end
        
        function declareGet(obj,TCsource,functionName,helpmsg)
        % declareGet(obj,TCsource,functionName) 
        %
        % Declares a 'get' operation to be compiled.  Each 'get'
        % operation retrives data from a compiled sparse
        % variables. TCsource may be a cell array for multiple
        % simultaneous gets.
            
            if ~iscell(TCsource) 
                TCsource={TCsource};
            end
            if isempty(TCsource)
                TCsource
                error('empty declareGet')
            end
            k=zeros(length(TCsource),1);
            for i=1:length(TCsource)
                k(i)=addTCexpression(obj,TCsource{i});
                %computeScalarInstructions(obj,k(i));
            end
            obj.gets(end+1,1)=...
                struct('functionName',{functionName},'source',{k},'parentGroups',{[]});
        
            % add to template
            obj.template(end+1,1).MEXfunction=sprintf('%s',functionName);
            obj.template(end,1).Sfunction=sprintf('%sS',functionName);
            obj.template(end,1).Cfunction=sprintf('%s',functionName);
            obj.template(end,1).method=sprintf('%s',functionName);
            if nargin<4
                obj.template(end,1).help='';
            else
                obj.template(end,1).help=helpmsg;
            end            
            for i=1:length(TCsource)
                obj.template(end,1).outputs(i).type=obj.scratchbookType;
                if strcmp(type(TCsource{i}),'variable')
                    obj.template(end,1).outputs(i).name=name(TCsource{i});
                else
                    obj.template(end,1).outputs(i).name=sprintf('output%d',i);
                end
                obj.template(end,1).outputs(i).sizes=TCsource{i}.osize;
            end
        end
        
        function declareCopy(obj,TCdestination,TCsource,functionName,helpmsg)
        % declareCopy(obj,TCdestination,TCsource,functionName)
        %
        % Declares a 'copy' operation to be compiled.  Each 'copy'
        % operation assigns data from one compiled sparse variable to
        % another one.  TCsource and TC destination may be cell
        % arrays for multiple simultaneous copies.
            
            if ~iscell(TCdestination) 
                TCdestination={TCdestination};
            end
            if ~iscell(TCsource) 
                TCsource={TCsource};
            end
            if length(TCdestination)~=length(TCsource)
                error('declareCopy: mismatch between lengths of source/destination cell arrays\n');
            end

            % remove 0-size variables
            for i=length(TCsource):-1:1
                if prod(size(TCsource{i}))==0
                    %fprintf('discarding copy of 0-sized variable \n');
                   TCsource(i)=[];
                   TCdestination(i)=[];
                end
            end            
            
            if isempty(TCdestination)
                TCsource
                TCdestination
                error('empty declareCopy')
            end
            
            % add expressions
            k1=zeros(length(TCdestination),1);
            k2=zeros(length(TCsource),1);            
            
            for i=1:length(TCsource)
                if ~strcmp(type(TCdestination{i}),'variable')
                    fprintf('declareCopy: ATTENTION destination ''%s'' is not a Tvariable, which may lead to unpredictable behavior\n',str(TCdestination{i}));
                end
                
                if ~isequal(size(TCdestination{i}),size(TCsource{i}))
                    disp('source =')
                    disp(TCsource{i})
                    disp('destination =');
                    disp(TCdestination{i})
                    error('declareCopy: source size [%s] differs from destination size [%s]\n',...
                          index2str(size(TCsource{i})),...
                          index2str(size(TCdestination{i})));
                end
                before=length(obj.vectorizedOperations);
                k1(i)=addTCexpression(obj,TCdestination{i});
                if length(obj.vectorizedOperations)>before
                    fprintf('declareCopy: ATTENTION setting a previously unseen variable %s, which may lead to unpredictable behavior\n',str(TCdestination{i}));
                end
                k2(i)=addTCexpression(obj,TCsource{i});
                %computeScalarInstructions(obj,k1(i));
                %computeScalarInstructions(obj,k2(i));
            end
            obj.copies(end+1,1)=struct('functionName',{functionName},'source',{k2},'destination',{k1},...
                                       'childrenGroups',{[]},'parentGroups',{[]});
        
        
            % add to template
            obj.template(end+1,1).MEXfunction=sprintf('%s',functionName);
            obj.template(end,1).Sfunction=sprintf('%sS',functionName);
            obj.template(end,1).Cfunction=sprintf('%s',functionName);
            obj.template(end,1).method=sprintf('%s',functionName);
            if nargin<5
                obj.template(end,1).help='';
            else
                obj.template(end,1).help=helpmsg;
            end            
        end
        
        % function [TCvar,subscripts,instructions]=declareComputation(obj,TCsource,name)
        % % [TCvar,subscripts,instructions]=declareComputation(obj,TCsource)
        % %
        % %   Adds a TC expression to the csparse obj, returning a
        % %   variable that serves as an alias to the expression and
        % %   also returning the subscripts of its nonzero elements.
        %     k=addTCexpression(obj,TCsource);
        %     %disp(obj)
        %     % if needed, compute instructions 
        %     subscripts=getOne(obj.vectorizedOperations,'subscripts',k);
        %     if isnan(subscripts)
        %         computeScalarInstructions(obj,k)
        %         subscripts=getOne(obj.vectorizedOperations,'subscripts',k);
        %     end
        %     instructions=getOne(obj.vectorizedOperations,'instructions',k);
        %     if nargin>=3 && ~isempty(name)
        %         % create linked variable
        %         TCvar=Tvariable(name,size(TCsource));
        %         old=height(obj.vectorizedOperations);
        %         k=addTCexpression(obj,TCvar);
        %         if height(obj.vectorizedOperations)==old
        %             error('declareComputation: linked variable ''%s'' already existed\n',name);
        %         end
        %         % link to existing instructions
        %         set(obj.vectorizedOperations,'subscripts',k,subscripts);
        %         set(obj.vectorizedOperations,'instructions',k,instructions);
        %     else
        %         TCvar=[];
        %     end
        %     %disp(obj)
        % end
        
        function [TCvar,subscripts,instructions]=declareAlias(obj,TCsource,name,atomic)
        % TCvar=declareAlias(obj,TCsource,name)
        % TCvar=declareAlias(obj,TCsource,name,atomic)
        % [TCvar,subscripts]=declareAlias(obj,TCsource,name,atomic)
        %
        %   Returns TC variable with the given name that acts as an
        %   alias for a TC expression, that is added to the csparse
        %   obj. 
        %   
        %   When atomic is true, the top expression is declared atomic
        %   
        %   Optionally, returns the TC variable's sparsity structure.
            
            if nargin<4
                atomic=false;
            end
            
            k=addTCexpression(obj,TCsource,atomic);
            
            if nargout>1
                subscripts=getOne(obj.vectorizedOperations,'subscripts',k);
                if isnan(subscripts)
                    computeScalarInstructions(obj,k)
                    subscripts=getOne(obj.vectorizedOperations,'subscripts',k);
                end
                instructions=getOne(obj.vectorizedOperations,'instructions',k);
            end
            
            if nargin>=3 && ~isempty(name)
                % create linked variable
                TCvar=Tvariable(name,size(TCsource));
                updateFile2table(TCvar,1);
                old=height(obj.vectorizedOperations);
                kV=addTCexpression(obj,TCvar,atomic);
                if height(obj.vectorizedOperations)==old
                    error('declareAlias: linked variable ''%s'' already existed\n',name);
                end
                % link to existing expression
                set(obj.vectorizedOperations,'operands',kV,k);
                if nargout>1
                    set(obj.vectorizedOperations,'subscripts',k,subscripts);
                    set(obj.vectorizedOperations,'instructions',k,instructions);
                end
            else
                TCvar=[];
            end
        end
        
        function declareSave(obj,TCsource,functionName,filename)
        % declareSave(obj,TCsource,functionName,filename)
        %   Writes to a file the subscripts of its nonzero elements, and
        %   declares a 'save' operation to be compiled. Each 'save'
        %   writes to a file the values of a compiled sparse variables
            k=addTCexpression(obj,TCsource);
            % magic=int64(intmax('int64')*rand(1)); 
            magic=int64(etime(clock(),[2000,1,1,0,0,0])*1e9); % ns 2000/1/1
            %fprintf('      declareSave(%s,magic=%d\n',functionName,magic);
            obj.saves(end+1,1)=struct('functionName',functionName,'filename',filename,...
                                      'source',k,'parentGroups',[],'magic',magic);
        end
        
        function declareFunction(obj,filename,functionName,defines,inputs,outputs,method,helpmsg)
        %   declarefunction(obj,filename,functionName,defines,inputs,outputs) 
        % or 
        %   declarefunction(obj,filename,functionName,defines) 
        %
        %   Declares a C (first form) or a matlab (second form) 
        %   function that typically calls the functions
        %   created through declareSet, declareGet, declareCopy,
        %   declareSave. The function can be found in the file
        %   'filename' and is called 'functionName()'. 
        %   (for matlab functions, the 1st function in the file
        %   will be changed to functionName(), if that was not the case)
        %
        %   For matlab functions, the structure 'defines' specifies a
        %   set of constants that will be included in the class and
        %   can be used to pass parameters to the matlab
        %   function, as in: 
        %     defines.name1 = value
        %     defines.name2 = value
        %
        %   For C functions, the inputs and outputs are passed by reference and
        %   are decribed in the structure arrays 'inputs' and
        %   'outputs', with fields
        %      .name = string with the name of the parameter
        %      .type = string with the matlab type of the
        %              parameter, as in
        %              (uint8|uint16|uint32|uint64|int8|int16|int32|int64|float|double)
        %      .size = array with the dimensions of the matrix
        %   The structure 'defines' specifies a set of #define
        %   pre-processor directives that should precede the C
        %   function definition and can be used to pass (hardcoded)
        %   parameters to the C function, as in:
        %      defines.name1 = {string or scalar}
        %      defines.name2 = {string or scalar}

            obj.externalFunctions(end+1).functionName=functionName;
            obj.externalFunctions(end).fileName=filename;
            obj.externalFunctions(end).defines=defines;
            obj.externalFunctions(end).inputs=inputs;
            obj.externalFunctions(end).outputs=outputs;
        
            % add to template
            obj.template(end+1,1).MEXfunction=sprintf('%s',functionName);
            obj.template(end,1).Sfunction=sprintf('%sS',functionName);
            obj.template(end,1).Cfunction=sprintf('%s',functionName);
            obj.template(end,1).method=method;
            obj.template(end,1).inputs=inputs;
            obj.template(end,1).outputs=outputs;
            if nargin<8
                obj.template(end,1).help='';
            else
                obj.template(end,1).help=helpmsg;
            end            
        
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create table with vectorized operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function k=addTCexpression(obj,TCobj,atomic)
        % k=addTCexpression(obj,TCobj)
        %   Parses a tenscalc expression into a sequence of
        %   ''vectorizedOperations'' and adds them to the csparse
        %   object.
        %
        %   When atomic is true, the top expression is declared atomic

            global TCsymbolicExpressions;

            if nargin<3
                atomic=false;
            end
            
            TCobj=toCalculus(TCobj);

            %% pre-processing
           
            %TCobj=precompute(TCobj);

            typ=type(TCobj);
            switch (typ)
              case 'traceinv'
                ops=operands(TCobj);
                op1=Tcalculus(ops(1));
                TCobj=trace(op1\Teye(size(op1)));
                typ=type(TCobj);
              case 'inv'
                ops=operands(TCobj);
                op1=Tcalculus(ops(1));
                TCobj=op1\Teye(size(op1));
                typ=type(TCobj);
              case 'interpolate'
                ops=operands(TCobj);
                TCobj=interpolate(Tcalculus(ops(1)),Tcalculus(ops(2)),...
                                  Tcalculus(ops(3)),Tcalculus(ops(4)),...
                                  TCobj.value,true);
                typ=type(TCobj);
              case 'Ginterpolate'
                ops=operands(TCobj);
                TCobj=Ginterpolate(Tcalculus(ops(1)),Tcalculus(ops(2)),...
                                   Tcalculus(ops(3)),Tcalculus(ops(4)),...
                                   TCobj.value,true);
                typ=type(TCobj);
              case 'Hinterpolate'
                ops=operands(TCobj);
                TCobj=Hinterpolate(Tcalculus(ops(1)),Tcalculus(ops(2)),...
                                   Tcalculus(ops(3)),Tcalculus(ops(4)),...
                                   TCobj.value,true);
                typ=type(TCobj);
            end
            if obj.tprod2matlab && strcmp(typ,'tprod') 
                TCobj=tprod_tprod2matlab(TCobj);
                typ=type(TCobj);
            end
            if length(obj.TCindex2CSvectorized)>=TCobj.TCindex && ...
                    obj.TCindex2CSvectorized(TCobj.TCindex)~=0
                % allreadey added TCobj
                k=obj.TCindex2CSvectorized(TCobj.TCindex);
                return
            end
                
            % add children
            ops=operands(TCobj);

            for i=1:length(ops)
                ops(i)=addTCexpression(obj,Tcalculus(ops(i)));
                optype=getOne(obj.vectorizedOperations,'type',ops(1));
                if strcmp(optype,'variable')
                    opop=getOne(obj.vectorizedOperations,'operands',ops(1));
                    if ~isempty(opop)
                        % alias
                        ops(i)=opop;
                    end
                end
            end
                
            osize=size(TCobj);
            description=file_line(TCobj);

            operandsMatch=true;
            switch typ
              case 'variable'
                oname=parameters(TCobj);
                nameMatch=true;
                operandsMatch=false;
                parametersMatch=false;
                pars=[];
              case 'constant'
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                pars=parameters(TCobj);
                parametersMatch=true;
              case {'plus','min2','max2'}
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                pars=cell2mat(op_parameters(TCobj));
                parametersMatch=true;
              case {'ones','eye','zeros','reshape'}
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                pars=osize;
                parametersMatch=true;
              case {'clp','subsref','min','max',...
                    'all','any','cat','sum','repmat','abs',...
                    'times','mtimes','norm2','norm1','norminf',...
                    'rdivide','ctranspose','diag','full'}
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                pars=parameters(TCobj);
                parametersMatch=true;
              case {'tprod','tprod_matlab'}
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                pars=op_parameters(TCobj);
                pars={parameters(TCobj),pars{:}}';
                parametersMatch=true;
                
              case 'compose'
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                pars=parameters(TCobj);
                pars=char(pars{1}); % derivatives do not need to match
                                    % (derivative of exp may be shorter because
                                    % truncated by gradient)
                parametersMatch=true;
                
              case {'logdet','traceinv'}
                pars=[];
                parametersMatch=false;
                optype=getOne(obj.vectorizedOperations,'type',ops(1));
                opatomic=getOne(obj.vectorizedOperations,'atomic',ops(1));
                switch optype
                  case {'ldl'}
                    typ=[typ,'_ldl'];
                  case {'lu','lu_sym'}
                    typ=[typ,'_lu'];
                  otherwise;
                    error('unexpected operand for logdet ''%s'' this operator can only be applied to a matrix that has been factorized using ''ldl'', ''lu'', or ''lu_sym''\n',optype);
                end
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                
              case 'mldivide'
                pars=[];
                parametersMatch=false;
                optype=getOne(obj.vectorizedOperations,'type',ops(1));
                opatomic=getOne(obj.vectorizedOperations,'atomic',ops(1));
                switch optype
                  case {'lu','lu_sym'}
                    if ~opatomic
                        %% The object creation should probably be done inside Tcalculus
                        % Since   L U x = b <=> L y = b & U x = y
                        % we expand
                        %    mldivide(LU,b) = mldivide_u(LU,mldivide_l1(LU,b))
                        % where 
                        %    mldivide_l1(LU,b) - 1) computes L by 
                        %                           . extracting the strictly lower-triangular 
                        %                             entries of LU
                        %                           . adding 1 to the diagonal
                        %                        2) applies any required row-permutation to b 
                        %                           (from the LU factorization)
                        %                        3) solves L y = b
                        %    mldivide_u(LU,y)  - 1) computes U by 
                        %                           . extracting the (non-strict) upper-triangular 
                        %                             entries of LU
                        %                        2) solves U x = y
                        %                        3) applies any required column-permutation to x
                        %                           (from the LU factorization)
                        typ='mldivide_l1';
                        description1=sprintf('mldivide_l1 for %s',description);
                        oname=sprintf('%s%d',typ,height(obj.vectorizedOperations)+1);
                        nameMatch=false;
                        k=appendRowUnique(obj.vectorizedOperations,typ,true,...
                                          description1,false,oname,nameMatch,...
                                          osize,true,NaN,false,...
                                          ops,true,pars,false,NaN,false,atomic);
                        obj.nAddedVectorizedOperations=obj.nAddedVectorizedOperations+1;
                        typ='mldivide_u';
                        description=sprintf('mldivide_u for %s',description);
                        oname=sprintf('%s%d',typ,height(obj.vectorizedOperations)+1);
                        ops(2)=k;
                        % k=appendRowUnique(obj.vectorizedOperations,'mldivide_u',true,...
                        %                   description1,false,oname,nameMatch,...
                        %                   osize,true,NaN,false,...
                        %                   ops,true,pars,false,NaN,false,atomic);
                        % obj.nAddedVectorizedOperations=obj.nAddedVectorizedOperations+1;
                    else
                        % Make sure b is full in   L U x = b
                        typ='full';
                        description1=sprintf('full for %s',description);
                        oname=sprintf('%s%d',typ,height(obj.vectorizedOperations)+1);
                        nameMatch=false;
                        k=appendRowUnique(obj.vectorizedOperations,typ,true,...
                                          description1,false,oname,nameMatch,...
                                          osize,true,NaN,false,...
                                          ops(2),true,pars,false,NaN,false,atomic);
                        obj.nAddedVectorizedOperations=obj.nAddedVectorizedOperations+2;
                        typ='mldivide';
                        oname=sprintf('%s%d',typ,height(obj.vectorizedOperations)+1);
                        ops(2)=k;
                    end
                  case 'ldl'
                    %% The object creation should really be done inside Tcalculus
                    % Since   L D U x = b <=> L y = b & D z = y & U x = z
                    % we expand
                    %    mldivide(LDU,b) = mldivide_u1(LDU,mldivide_d(LDU,mldivide_l1(LDU,b)))
                    % where 
                    %    mldivide_l1(LDU,b) - 1) computes L by 
                    %                            . extracting the strictly lower-triangular 
                    %                              entries of LDU
                    %                            . adding 1 to the diagonal
                    %                         2) applies any required permutation to b 
                    %                            (from the LDU factorization)
                    %                         3) solves L y = b
                    %    
                    %    mldivide_d(LDU,y)  - 1) computes D by 
                    %                            . extracting the main diagonal of LDU
                    %                         2) solves D z = y
                    %    
                    %    mldivide_u1(LDU,z) - 1) computes U by 
                    %                            . extracting the strictly upper-triangular 
                    %                              entries of LU
                    %                            . adding 1 to the diagonal
                    %                         2) solves U x = z
                    %                         3) applies any required permutation to x
                    %                            (from the LDU factorization)
                    typ='mldivide_l1';
                    description1=sprintf('mldivide_l1 for %s',description);
                    oname=sprintf('%s%d',typ,height(obj.vectorizedOperations)+1);
                    nameMatch=false;
                    k=appendRowUnique(obj.vectorizedOperations,typ,true,...
                                      description1,false,oname,nameMatch,...
                                      osize,true,NaN,false,...
                                      ops,true,pars,false,NaN,false,atomic,true);
                    obj.nAddedVectorizedOperations=obj.nAddedVectorizedOperations+1;
                    typ='mldivide_d';
                    description1=sprintf('mldivide_d for %s',description);
                    oname=sprintf('%s%d',typ,height(obj.vectorizedOperations)+1);
                    ops(2)=k;
                    k=appendRowUnique(obj.vectorizedOperations,typ,true,...
                                      description1,false,oname,nameMatch,...
                                      osize,true,NaN,false,...
                                      ops,true,pars,false,NaN,false,atomic,true);
                    obj.nAddedVectorizedOperations=obj.nAddedVectorizedOperations+1;
                    typ='mldivide_u1';
                    description=sprintf('mldivide_u1 for %s',description);
                    oname=sprintf('%s%d',typ,height(obj.vectorizedOperations)+1);
                    ops(2)=k;
                    % k=appendRowUnique(obj.vectorizedOperations,'mldivide_u1',true,...
                    %                   description1,false,oname,nameMatch,...
                    %                   osize,true,NaN,false,...
                    %                   ops,true,pars,false,NaN,false,atomic,true);
                    % obj.nAddedVectorizedOperations=obj.nAddedVectorizedOperations+1;
                otherwise;
                    error('unexpected operand for mldivide ''%s'' this operator can only be applied to a matrix that has been factorized using ''ldl'', ''lu'', or ''lu_sym''\n',optype);
                end
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;

              case 'lu'
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                parametersMatch=false; % since matrix may already have p & q defined
                pars=parameters(TCobj);
                pars={[],[],pars.typical_subscripts,pars.typical_values};
                  
              case 'lu_sym'
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                parametersMatch=false; % since matrix may already have p & q defined
                pars=parameters(TCobj);
                pars={[],[],pars.typical_subscripts,pars.typical_values};
                  
              case 'ldl'
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                parametersMatch=false; % since matrix may already have p defined
                pars=parameters(TCobj);
                pars={[],[],pars.typical_subscripts,pars.typical_values};
                  
              case 'ldl_d'
                optype=getOne(obj.vectorizedOperations,'type',ops(1));
                if ~strcmp(optype,'ldl')
                    error('unexpected operand for ldl_d ''%s'' this operator can only be applied to a matrix that has been factorized using ''ldl''\n',optype);
                end                    
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                parametersMatch=false;
                pars=parameters(TCobj);
                pars=[];
                  
              case 'ldl_l'
                optype=getOne(obj.vectorizedOperations,'type',ops(1));
                if ~strcmp(optype,'ldl')
                    error('unexpected operand for ldl_l ''%s'' this operator can only be applied to a matrix that has been factorized using ''ldl''\n',optype);
                end                    
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                parametersMatch=false;
                pars=parameters(TCobj);
                pars=[];
                  
              case 'chol'
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                parametersMatch=false; % since matrix may already have p defined
                pars=parameters(TCobj);
                pars={[],[],pars.typical_subscripts,pars.typical_values};
              
              otherwise
                error('addTCexpression: unimplemented type %s\n',typ)
                oname=sprintf('%s_%d',char(typ),height(obj.vectorizedOperations)+1);
                nameMatch=false;
                parametersMatch=false;
                pars=[];
            end
            
            % add this expression
            len=height(obj.vectorizedOperations);

            k=appendRowUnique(obj.vectorizedOperations,typ,true,...
                              description,false,oname,nameMatch,...
                              osize,true,NaN,false,...
                              ops,operandsMatch,pars,parametersMatch,NaN,false,atomic,true);
            obj.TCindex2CSvectorized(TCobj.TCindex)=k;
            obj.nAddedVectorizedOperations=obj.nAddedVectorizedOperations+1;

            if mod(obj.nAddedVectorizedOperations,10000)==0
                fprintf('added %d TC expr (%d unique)...',obj.nAddedVectorizedOperations,len);
            end

            %fprintf('addTCexpression: %-20s stored at row %4d, nAddedVectorizedOperations=%6d, table height=%4d\n',['"',oname,'"'],k,obj.nAddedVectorizedOperations,len);
            
            if obj.debug>=3 && k>len
                if ~strcmp(type(TCobj),'full');
                    addTCexpression(obj,full(TCobj));
                else
                    oname=sprintf('debugGet_%s',...
                                  getOne(obj.vectorizedOperations,'name',ops(1)));
                    obj.gets(end+1,1)=struct(...
                        'functionName',{name},'source',{k},'parentGroups',{[]});
                end
            end
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create table with instruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function instr=newInstructions(obj,types,parameters,operands,vectorizedOperation,fast)
        % Creates one or several instructions of the given type, with
        % the given parameters, operands's intructions, and
        % corresponding vectorizedOperation.  If a similar intruction
        % already exists, returns the previous location.
        %
        % The input parameter 'type' may be a string or a cell array.
        % The input parameters 'parameters' and 'operands', must be
        % cell arrays (one per instruction to be created). 
            if nargin<6
                fast=false;  % appendUnique needed to detect structural symmetries
            end
            if isnumeric(types)
                types={types};
            end
            n=max([length(types),length(parameters),length(operands)]);
            if length(types)<n
                types=repmat(types,n,1);
            end
            if length(parameters)<n
                parameters=repmat(parameters,n,1);
            end
            if length(operands)<n
                operands=repmat(operands,n,1);
            end
            instr=nan(n,1);
            for i=1:n
                % not clear why exception for I_set (2015/03/05)
                if types{i}==obj.Itypes.I_set || fast  
                    instr(i)=appendInstruction(types{i},parameters{i},int64(operands{i}));
                else
                    [operands{i},parameters{i}]=normalize(obj,types{i},parameters{i},operands{i});
                    instr(i)=appendUniqueInstruction(types{i},parameters{i},int64(operands{i}));
                end
                obj.statistics.nAddedInstructions=obj.statistics.nAddedInstructions+1;
                
                if mod(obj.statistics.nAddedInstructions,100000)==0
                    fprintf('added %d instr. (%d unique)...',...
                            obj.statistics.nAddedInstructions,instructionsTableHeight());
                end
            end

        end

        function instr=newInstruction(obj,typ,parameters,operands,vectorizedOperation,fast)
        % Creates a single instruction of the given type, with the given
        % parameters, operands, and corresponding vectorizedOperation.
        % If a similar intruction already exists, returns the
        % previous location. 
        %
        % Slightly faster than newInstructions for a single instruction.
            if nargin<6
                fast=false; % appendUnique needed to detect structural symmetries
            end
            % not clear why exception for I_set (2015/03/05)
            if typ==obj.Itypes.I_set || fast
                instr=double(appendInstruction(typ,parameters,int64(operands)));
            else
                [operands,parameters]=normalize(obj,typ,parameters,operands);
                instr=double(appendUniqueInstruction(typ,parameters,int64(operands)));
            end
            obj.statistics.nAddedInstructions=obj.statistics.nAddedInstructions+1;
            
            if mod(obj.statistics.nAddedInstructions,100000)==0
                 fprintf('added %d instr. (%d unique)...',...
                         obj.statistics.nAddedInstructions,instructionsTableHeight());
            end
        end

        function [operands,parameters]=normalize(obj,typ,parameters,operands)
            switch(typ)
              case obj.Itypes.I_sumprod
                if ~isequal(size(operands),parameters)
                    operands=reshape(operands,parameters);
                end
                operands=sort(operands,1);     % sort factors in each product;
                operands=sortrows(operands')'; % sort sums 
              case obj.Itypes.I_sum
                [operands,k]=sort(operands);
                parameters=parameters(k);
              case obj.Itypes.I_minus_dot
                if length(operands)>1
                    operands=reshape(operands,2,[]);
                    operands=sort(operands,1);     % sort factors in each product;
                    operands=sortrows(operands')'; % sort sums 
                    operands=reshape(operands,1,[]);
                end
              case obj.Itypes.I_minus_dot_div
                if length(operands)>2
                    ops=reshape(operands(2:end),2,[]);
                    ops=sort(ops,1);     % sort factors in each product;
                    ops=sortrows(ops')'; % sort sums 
                    operands(2:end)=reshape(ops,1,[]);
                end
              case obj.Itypes.I_plus_minus_dot
                if length(operands)>2
                    ops=reshape(operands(2:end),2,[]);
                    ops=sort(ops,1);     % sort factors in each product;
                    ops=sortrows(ops')'; % sort sums 
                    operands(2:end)=reshape(ops,1,[]);
                end
              case obj.Itypes.I_plus_minus_dot_div
                if length(operands)>3
                    ops=reshape(operands(3:end),2,[]);
                    ops=sort(ops,1);     % sort factors in each product;
                    ops=sortrows(ops')'; % sort sums 
                    operands(3:end)=reshape(ops,1,[]);
                end
              case obj.Itypes.I_plus_sqr
                operands=sort(operands);
              case obj.Itypes.I_min
                operands=sort(operands);
              case obj.Itypes.I_min0
                operands=sort(operands);
              case obj.Itypes.I_max
                operands=sort(operands);
              case obj.Itypes.I_max0
                operands=sort(operands);
              case obj.Itypes.I_max_abs
                operands=sort(operands);
              case obj.Itypes.I_clp
                if length(operands)>2
                    operands=reshape(operands,2,[]);
                    operands=sortrows(operands')'; % sort inequalities
                    operands=reshape(operands,1,[]);
                end
              otherwise
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            function children=childrenOf(obj,parents)
            % Computes the set of instructions that depend on 'parents' (exclusive)
            
            children=obj.dependencyGraph*sparse(parents,ones(size(parents)),ones(size(parents)),size(obj.dependencyGraph,1),1);
            while 1
                new=double((children+obj.dependencyGraph*children)>0);
                %[(1:length(children))',full(children),full(new)]
                if isequal(new,children)
                    break
                end
                children=new;
            end
        end

        function parents=parentsOf(obj,children)
            % Computes the set of instructions that the 'children' depend on (inclusive)
            parents=sparse(ones(size(children)),children,ones(size(children)),1,size(obj.dependencyGraph,1));
            while 1
                new=double((parents+parents*obj.dependencyGraph)>0);
                if isequal(new,parents)
                    break
                end
                parents=new;
            end
            parents=parents';
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compile code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function compile2C(obj,minInstructions4loop,maxInstructionsPerFunction,Cfunction,Hfunction,logFile,folder,profiling)
        % compile2C(obj,minInstructions4loop,maxInstructionsPerFunction,Cfunction,Hfunction,logFile,folder,profiling)
        %
        % Compiles code to C and appends it to a given file.
        %
        % Inputs:
        % obj - csparse object
        % minInstructions4loop - minimum number of similar instructions to
        %                        be implmented as a for loop (rather than inlined) 
        % maxInstructionsPerFunction - maximum number of instructions to be
        %                        included in a single function. Used to avoid
        %                        very large functions for which the compiler
        %                        could misbehave.
        % Cfunction - file where the C code should be written
        % Hfunction - file where the C function headers should be written (optional)
        % logFile   - file where statistics information should be written (optional)
        % folder    - folder where all files will be written (optional)
        % profiling - when non-zero adds profiling code (optional)
            
            if nargin<8
                profiling=false;
            end
            if nargin<7
                folder='.';
            end
            if nargin<6
                logFile='';
            end
            if nargin<5
                Hfunction='';
            end

            fprintf('  computeScalarInstructions...\n');t0=clock;
            checkVariableSets(obj);
            computeScalarInstructions(obj,1:height(obj.vectorizedOperations),folder);
            fprintf('    done (%.3f sec)\n',etime(clock(),t0));
            fprintf('  dependencyGroups... ');t0=clock;
            dependencyGroups(obj);
            fprintf('done computeScalarInstructions (%.3f sec)\n',etime(clock(),t0));

            fprintf('  compile2C: %d added vectorized Operations, %d unique vectorized Operations\n',...
                    obj.nAddedVectorizedOperations,height(obj.vectorizedOperations));
            fprintf('  compile2C: %d added instructions, %d unique instructions\n',...
                    obj.statistics.nAddedInstructions,instructionsTableHeight());
            fprintf('  compile2C: sizeScratchbook=%d %s, nGroups=%d, nSets=%d, nGets=%d, nCopies=%d\n',...
                    max(obj.memoryLocations),obj.scratchbookType,...
                    size(obj.dependencyGroups,1),...
                    length(obj.sets),length(obj.gets),length(obj.copies));

            fprintf('  write C code... ');t0=clock;
            %writeCswitchpergroup(obj,codeType,Cfunction,Hfunction,logFile,profiling);
            writeCfunctionpergroup(obj,minInstructions4loop,maxInstructionsPerFunction,...
                                   Cfunction,Hfunction,logFile,folder,profiling);
            fprintf('done (%.3f sec)\n',etime(clock(),t0));
        end
    
        function compile2matlab(obj,Mfunction,logFile,classhelp,profiling)
        % compile2matlab(obj,Mfunction,logFile,classhelp,profiling)
        %
        % Compiles code to Matlab and appends it to a given file
        % Inputs:
        % obj - csparse object
        % Mfunction - file where the M code should be written
        % logFile   - file where statistics information should be written (optional)
        % profiling - when non-zero adds profiling code (optional)
            if nargin<4
                profiling=false;
            end
            if nargin<3
                logFile='';
            end

            fprintf('  computeMatlabInstructions... ');t0=clock;
            checkVariableSets(obj);
            computeMatlabInstructions(obj);
            fprintf('done (%.3f sec)\n',etime(clock(),t0));
            fprintf('  dependencyGroups... ');t0=clock;
            dependencyGroups(obj);
            fprintf('done (%.3f sec)\n',etime(clock(),t0));

            fprintf('  compile2matlab: %d added vectorized Operations, %d unique vectorized Operations\n',...
                    obj.nAddedVectorizedOperations,height(obj.vectorizedOperations));
            fprintf('  compile2matlab: %d added instructions, %d unique instructions\n',...
                    obj.statistics.nAddedInstructions,instructionsTableHeight());
            fprintf('  compile2matlab: sizeScratchbook=%d %s, nGroups=%d, nSets=%d, nGets=%d, nCopies=%d\n',...
                    max(obj.memoryLocations),obj.scratchbookType,...
                    size(obj.dependencyGroups,1),...
                    length(obj.sets),length(obj.gets),length(obj.copies));

            writeMatlab(obj,Mfunction,logFile,classhelp);
            fprintf('done (%.3f sec)\n',etime(clock(),t0));
        end
        
        function checkVariableSets(obj)
        % checkVariableSets(obj)
        %
        % Makes sure that every variable in the csparse object has
        % a corresponding set or copy
            
            variables=[];
            variable_names={};
            operands=[];
            for thisExp=1:height(obj.vectorizedOperations)
                typ=getOne(obj.vectorizedOperations,'type',thisExp);
                ops=getOne(obj.vectorizedOperations,'operands',thisExp);
                if strcmp(typ,'variable') && isempty(ops) % variable & not alias
                    variables(end+1,1)=thisExp;
                    name=getOne(obj.vectorizedOperations,'name',thisExp);
                    variable_names{end+1,1}=name;
                    %fprintf('%d: type="%s", name="%s", operands=%d\n',thisExp,typ,name,ops);
                end
            end
            variables=unique(variables);
            
            sets=[];
            set_names={};
            for k=1:length(obj.sets)
                sets=[sets;obj.sets(k).destination];
                set_names{end+1,1}=obj.sets(k).functionName;
            end
            for k=1:length(obj.copies)
                sets=[sets;obj.copies(k).destination];
                set_names{end+1,1}=obj.copies(k).functionName;
            end
            sets=unique(sets);
            
            missing=setdiff(variables,sets);
            if ~isempty(missing)
                fprintf('\n\n');
                for i=1:length(missing)
                    thisExp=missing(i);
                    typ=getOne(obj.vectorizedOperations,'type',thisExp);
                    variables(end+1,1)=thisExp;
                    name=getOne(obj.vectorizedOperations,'name',thisExp);
                    fprintf('%d: type="%s", name="%s"\n',thisExp,typ,name);
                end
                    
                error('variables missing sets/copies');
            end
            if ~isequal(variables,sets)
                variables
                variable_names
                sets
                set_names
                error('variables do not matchs sets/copies');
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Debug stuff
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function showInstruction(obj,inst,indent)
            [typ,parameters,operands]=getInstruction(int64(inst));
            if length(operands)>0
                fn=fields(obj.Itypes);
                fprintf([repmat(' ',1,indent),'%d %s (%d): '],inst,fn{typ},typ);
                fprintf('%d ',operands);
                fprintf('\n');
                for i=1:length(operands)
                    showInstruction(obj,operands(i),indent+1)
                end
            end
        end

        function [subs,instr]=getSparsity(obj,TCobj)
            k=addTCexpression(obj,TCobj);
            computeScalarInstructions(obj);
            %disp(obj)
            subs=getOne(obj.vectorizedOperations,'subscripts',k);
            if nargout>1
                instr=cell(size(subs,2),1);
                ks=getOne(obj.vectorizedOperations,'instructions',k);
                for k=1:length(ks)
                    [instr{k},parameters,o]=getInstructions(obj.instructions,int32(ks(k)));
                    if instr{k}==obj.Itypes.I_load;
                        instr{k}=sprintf('%s %g\n',instr{k},parameters);
                    end
                end
            end
        end
    
    end
    
end 

