function writeMatlab(obj,Mfunction,logFile,classHelp)
% writeMatlab(obj,Mfunction)
%    Write matlab code to implements all the gets, sets, and copies in
%    a scparse object.
% Inputs:
% obj - csparse object
% Mfunction - file where the Matlab code should be written
% classHelp - this string or string array is included in the
%             classdef file to provide an help message.
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

verboseLevel=0;

if verboseLevel>0
    fprintf('  writeMatlab ');
    t0=clock;
end

if ~isempty(logFile)
    fil=fopen(logFile,'w');
end

[~,className,~]=fileparts(Mfunction);

computeFunction='compute';
callComputeFunction=sprintf('      %s_%%d(obj);\n',computeFunction);

fid=fopen(Mfunction,'w');
%groups=getMulti(obj.instructions,'group');
groups=obj.instructionsGroup;

nGroups=max(groups);

if verboseLevel>0
    fprintf('WARNING: writeMatlab: each instruction is being assigned a specific memory address\n\tbut one should reuse memory\n');
end

fprintf(fid,'classdef %s < handle\n',className);
if ischar(classHelp)
    classHelp{1}=classHelp;
end
for i=1:length(classHelp)
    fprintf(fid,'%% %s\n',classHelp{i});
end
fprintf(fid,'  properties\n');

%% Parameters
fprintf(fid,'    %% defines\n');
for i=1:length(obj.externalFunctions)
    names=fields(obj.externalFunctions(i).defines);
    for j=1:length(names)
        value=getfield(obj.externalFunctions(i).defines,names{j});
        if ischar(value)
            fprintf(fid,'    %s=''%s'';\n',names{j},value);
        elseif isnumeric(value) 
            fprintf(fid,'    %s=[%s];\n',names{j},mymat2str(value));
        else
            value
            error('define %s has invalid value\n',names{j});
        end
    end
end

%% Storage areas
fprintf(fid,'    %% scrapbook\n');
fprintf(fid,'    m%d=[];\n',unique(obj.memoryLocations));
fprintf(fid,'    %% status of dependency groups\n');
fprintf(fid,'    groupStatus=[%s];\n',index2str(zeros(1,nGroups)));
fprintf(fid,'  end %% properties\n');
fprintf(fid,'  methods\n');

if ~isempty(logFile)

    if 1
        %% Write description of the csparse object
        fprintf(fil,'* Sparse object:\n');
        fprintf(fil,'%s',str(obj,1));
    end

    % c=[obj.TCindex2CSvectorized,(1:length(obj.TCindex2CSvectorized))'];
    % c(obj.TCindex2CSvectorized==0,:)=[];
    % [~,k]=sort(c(:,1));
    % c=c(k,:);
    % fprintf(fil,'CSparse index, TC index\n');
    % fprintf(fil,'  %6d %6d\n',c');
    
    %% Write Dependency groups to the log file
    fprintf(fil,'\n* Dependency groups:\n');
    fprintf(fil,'%35s','group');
    fprintf(fil,'%5d',0:nGroups-1);
    fprintf(fil,'\n');
    fprintf(fil,'%s\n',repmat('-',1,35));
    for i=1:size(obj.dependencyGroups,2)
        fprintf(fil,'%-35s',obj.dependencyGroupColName{i});
        fprintf(fil,'%5d',full(obj.dependencyGroups(:,i)));
        fprintf(fil,'\n');
    end
    fprintf(fil,'\n');
end

if verboseLevel>0
    fprintf('writeMatlab: creating computation functions\n');
end

%% Create functions to do the computations
for i=1:nGroups
    functionName=sprintf('%s_%d',computeFunction,i);
    if length(functionName)>namelengthmax
        error('function name is too long "%s"',functionName);
    end
    if verboseLevel>0
        fprintf('    function %s(obj)\n',functionName);
    end
    fprintf(fid,'    function %s(obj)\n',functionName);
    k=find(groups==i);
    fprintf(fid,'       if (obj.groupStatus(%d)==0)\n',i);
    if obj.debug>0
        fprintf(fid,'         fprintf(''computing group %d\\n'');\n',i);
    end
    writeMatlabInstructions(obj,fid,k);
    fprintf(fid,'         obj.groupStatus(%d)=1;\n',i);
    if obj.debug>0
        fprintf(fid,'       else fprintf(''(skipping group %d)\\n'');\n',i);
    end
    fprintf(fid,'       end\n');
    fprintf(fid,'    end\n');
end

if verboseLevel>0
    fprintf('writeMatlab: creating methods functions\n');
end

%% Create functions for sets
if verboseLevel>0
    fprintf('writeMatlab: creating sets functions\n');
end
for i=1:length(obj.sets)
    if length(obj.sets(i).functionName)>namelengthmax
        error('function name is too long "%s"',obj.sets(i).functionName);
    end
    if verboseLevel>0
        fprintf('    function %s(obj,input);\n',obj.sets(i).functionName);
    end
    fprintf(fid,'    function %s(obj,input)\n',obj.sets(i).functionName);
    if obj.debug>0
        fprintf(fid,'      fprintf(''running %s\\n'');\n',obj.sets(i).functionName);
    end
    instructions=getOne(obj.vectorizedOperations,'instructions',obj.sets(i).destination);
    osize=getOne(obj.vectorizedOperations,'osize',obj.sets(i).destination);
    while length(osize)<2
        osize(end+1)=1;
    end
    fprintf(fid,'      if ~isequal(size(input),[%s]),error(''%s: size mismatch: expected=[%s] received=[%%s]\\n'',index2str(size(input))),end;\n',...
            index2str(osize),obj.sets(i).functionName,index2str(osize));
    fprintf(fid,'      obj.m%d=input;\n',obj.memoryLocations(instructions));
    % set group status
    if isempty(obj.sets(i).childrenGroups)
        %fprintf('writeMatlab: %s has no effect\n',obj.sets(i).functionName);
    else
        fprintf(fid,'      obj.groupStatus(%d)=0;\n',obj.sets(i).childrenGroups);
    end
    fprintf(fid,'    end\n');
end

%% Create functions for gets
if verboseLevel>0
    fprintf('writeMatlab: creating gets functions\n');
end
for i=1:length(obj.gets)
    nSources=length(obj.gets(i).source);
    tmpl=obj.template(obj.gets(i).templateNdx);
    % write function header
    if verboseLevel>0
        fprintf('    function ');
    end
    fprintf(fid,'    function ');
    sep='[';
    for j=1:nSources
        if verboseLevel>0
            fprintf('%c%s',sep,tmpl.outputs(j).name);
        end
        fprintf(fid,'%c%s',sep,tmpl.outputs(j).name);
        sep=',';
    end
    if length(obj.gets(i).functionName)+1>namelengthmax
        error('function name is too long "%s"',obj.gets(i).functionName);
    end
    if verboseLevel>0
        fprintf(']=%s(obj)\n',obj.gets(i).functionName);
    end
    if nSources>1
        fprintf(fid,']=%s_(obj)\n',obj.gets(i).functionName);
    else
        fprintf(fid,']=%s(obj)\n',obj.gets(i).functionName);
    end
    % write function body
    if obj.debug>0
        fprintf(fid,'      fprintf(''running %s\\n'');\n',obj.gets(i).functionName);
    end
    
    % perform needed computations
    if ~isempty(obj.gets(i).parentGroups)
        fprintf(fid,callComputeFunction,obj.gets(i).parentGroups);
    end
    
    for j=1:nSources
        instructions=getOne(obj.vectorizedOperations,'instructions',obj.gets(i).source(j));
        fprintf(fid,'      %s=obj.m%d;\n',tmpl.outputs(j).name,obj.memoryLocations(instructions));
    end
    fprintf(fid,'    end\n');
    if nSources>1
        % output as structure
        if length(obj.gets(i).functionName)+7>namelengthmax
            error('function name is too long "%s"',obj.gets(i).functionName);
        end
        fprintf(fid,'    function varargout=%s(obj)\n',obj.gets(i).functionName);
        fprintf(fid,'       if nargout>1\n');
        fprintf(fid,'         varargout=cell(nargout,1);[varargout{:}]=%s_(obj);\n',obj.gets(i).functionName);
        fprintf(fid,'       else\n');
        fprintf(fid,'         ');
        sep='[';
        for j=1:nSources
            fprintf(fid,'%cvarargout{1}.%s',sep,tmpl.outputs(j).name);
            sep=',';
        end
        fprintf(fid,']=%s_(obj);\n',obj.gets(i).functionName);
        fprintf(fid,'       end\n');
        fprintf(fid,'    end\n');
    end
end

%% Create functions for copies
if verboseLevel>0
    fprintf('writeMatlab: creating copies functions\n');
end
for i=1:length(obj.copies)
    if length(obj.copies(i).functionName)>namelengthmax
        error('function name is too long "%s"',obj.copies(i).functionName);
    end
    % write function header
    if verboseLevel>0
        fprintf('    function %s(obj)\n',obj.copies(i).functionName);
    end
    fprintf(fid,'    function %s(obj)\n',obj.copies(i).functionName);
    if obj.debug>0
        fprintf(fid,'      fprintf(''running %s\\n'');\n',obj.copies(i).functionName);
    end
    
    % perform needed computations
    if ~isempty(obj.copies(i).parentGroups)
        fprintf(fid,callComputeFunction,obj.copies(i).parentGroups);
    end
    
    for k=1:length(obj.copies(i).destination)
        instructionsDestination=getOne(obj.vectorizedOperations,...
                                       'instructions',obj.copies(i).destination(k));
        instructionsSource=getOne(obj.vectorizedOperations,...
                                  'instructions',obj.copies(i).source(k));
        
        fprintf(fid,'      obj.m%d=obj.m%d;\n',...
                obj.memoryLocations(instructionsDestination),...
                obj.memoryLocations(instructionsSource));
    end
    
    if ~isempty(obj.copies(i).childrenGroups)
        % set group status
        fprintf(fid,'      obj.groupStatus(%d)=0;\n',obj.copies(i).childrenGroups);
    end
    
    fprintf(fid,'    end\n');
end

%% Create functions for externalFunctions
if verboseLevel>0
    fprintf('writeMatlab: creating external functions\n');
end
for i=1:length(obj.externalFunctions)
    if length(obj.externalFunctions(i).functionName)>namelengthmax
        error('function name is too long "%s"',obj.externalFunctions(i).functionName);
    end
    if verboseLevel>0
        fprintf('    function %s(obj,...)\n',obj.externalFunctions(i).functionName);
    end
    fprintf(fid,'    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'    %% Start of code for function %s(), copied from \"%s\" */\n',...
            obj.externalFunctions(i).functionName,obj.externalFunctions(i).fileName);
    fprintf(fid,'    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

    % append function
    fic=fopen(obj.externalFunctions(i).fileName,'r');
    rep=true;    search='function\s*([^=]*=\s*|)(\w+)(';
    replac=sprintf('function $1%s(',obj.externalFunctions(i).functionName);
    while 1
        tline=fgetl(fic);
        if ~ischar(tline)
            break
        end
        if rep && ~isempty(regexp(tline,search))
            tline=regexprep(tline,search,replac);
            rep=false;
        end
        fprintf(fid,'%s\n',tline);
    end
    fclose(fic);
    fprintf(fid,'\n');

    
    fprintf(fid,'    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'    %% End of code for function %s(), copied from \"%s\"\n',...
            obj.externalFunctions(i).functionName,obj.externalFunctions(i).fileName);
    fprintf(fid,'    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

end

%% Close

fprintf(fid,'  end %% methods\n');
fprintf(fid,'end %% classdef\n\n');
includeFile(fid,'GPL.m');
fprintf(fid,'%%/***********************************************/\n');
fprintf(fid,'%%/*  End of code generated by compile2Matlab.m  */\n');
fprintf(fid,'%%/***********************************************/\n\n');

fclose(fid);

rehash path;

end

function includeFile(fid,filename)
% includeFile(fid,filename)
%   includes the file named 'filename' into the file with the given
%   descriptor
    
    fii=fopen(filename,'r');
    if fii<0
        error('unable to open file ''%s''\n',filename);
    end
    fprintf(fid,'%% START OF #included "%s"\n',filename);
    a=fread(fii,inf);
    fwrite(fid,a);
    fprintf(fid,'%% END OF #included "%s"\n',filename);
    fclose(fii);
    
end
