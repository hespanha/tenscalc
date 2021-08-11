function dependencyGroups(obj)
% dependencyGroups(obj)
%
%   Clusters the instructions into "dependency groups" so that the
%   sets of
%     1) all the instructions that need to be recomputed due to a set
%        (or copy)
%     2) all the instructions that need to be computed to generate the
%        variable in a get (or copy)
%   can be expressed as unions of "dependency groups
%
%   This facilitates determining (in run time) which instructions need
%   to be executed in response to gets/sets/copies since one only
%   needs to keep track of each dependency groups need to be executed,
%   rather than which individual instructions need to be recomputed.
%
%   Attention: currently no attention is paid to the order in which
%   instructions (or dependency groups) need to be computed.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

verboseLevel=0;

findUnusedMemory=true;
reuseMemory=true;
compressMemory=true;
reuseLocalMemory=false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute dependency graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nInstructions=double(instructionsTableHeight());

%% Get dependencies from csparse/instructionsTable.c cmex functions
[children,parents]=getDependencies();
obj.dependencyGraph=sparse(double(children),double(parents),...
                           ones(size(children)),...
                           nInstructions,nInstructions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute dependency groups from gets/sets/copies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computes dependencies for gets, sets, copies

% create tranpose matrix since it is faster to assign rows of sparse matrices rather than columns

t0=clock;
fprintf(' dependency matrix with %d instructions',nInstructions);
dependencies=sparse([],[],[],...
                    length(obj.gets)+length(obj.sets)+2*length(obj.copies),nInstructions);

inputInstructions=zeros(1,0);
outputInstructions=zeros(1,0);
obj.dependencyGroupColName={};
j=1;
fprintf('(gets');
for i=1:length(obj.gets)
    instr=getMulti(obj.vectorizedOperations,'instructions',obj.gets(i).source);
    instr=vertcat(instr{:});
    dependencies(j,:)=parentsOf(obj,instr);
    obj.dependencyGroupColName{j}=obj.gets(i).functionName;
    j=j+1;
    outputInstructions=union(outputInstructions,instr);
end
fprintf(',saves');
for i=1:length(obj.saves)
    instr=getMulti(obj.vectorizedOperations,'instructions',obj.saves(i).source);
    instr=vertcat(instr{:});
    dependencies(j,:)=parentsOf(obj,instr);
    obj.dependencyGroupColName{j}=obj.saves(i).functionName;
    j=j+1;
    outputInstructions=union(outputInstructions,instr);
end
fprintf(',copies');
for i=1:length(obj.copies)
    instr=getMulti(obj.vectorizedOperations,'instructions',obj.copies(i).source);
    instr=vertcat(instr{:});
    dependencies(j,:)=parentsOf(obj,instr);
    obj.dependencyGroupColName{j}=obj.copies(i).functionName;
    j=j+1;
    outputInstructions=union(outputInstructions,instr);
end
fprintf(',sets');
for i=1:length(obj.sets)
    instr=getMulti(obj.vectorizedOperations,'instructions',obj.sets(i).destination);
    instr=vertcat(instr{:});
    dependencies(j,:)=childrenOf(obj,instr);
    obj.dependencyGroupColName{j}=obj.sets(i).functionName;
    j=j+1;
    inputInstructions=union(inputInstructions,instr);
end
fprintf(',copies');
for i=1:length(obj.copies)
    instr=getMulti(obj.vectorizedOperations,'instructions',obj.copies(i).destination);
    instr=vertcat(instr{:});
    dependencies(j,:)=childrenOf(obj,instr);
    obj.dependencyGroupColName{j}=obj.copies(i).functionName;
    j=j+1;
    inputInstructions=union(inputInstructions,instr);
end

% back to non-transpose matrix
dependencies=dependencies';
fprintf(' %.2f sec) determining groups',etime(clock,t0));

%% Determine unique dependency groups
groups=zeros(nInstructions,1);
% I_set instructions do not perform any action, so should be ignored
k=~findInstructionsByType(obj.Itypes.I_set);
[dependencyGroups,~,kGroup]=unique(dependencies(k,:),'rows');
groups(k)=kGroup;

%full(dependencies(k,:))
%full(dependencyGroups)

nGroups=size(dependencyGroups,1);

%% Break non-contiguous groups -- this should make sure there are no
%% cyclic dependencies among dependency groups
%% However, it may create more groups than are really needed, leading to extra overhead.
if 0
    %[(1:nInstructions)',groups],oldgroups=groups;
    for g=1:nGroups
        k=find(groups==g);
        br=find(k(2:end)~=k(1:end-1)+1);
        if ~isempty(br)
            if verboseLevel>0
                fprintf('  dependencyGroup: break group %d in %d parts\n',g,length(br)+1);
            end
            br(end+1)=length(k);
            for i=1:length(br)-1
                nGroups=nGroups+1;
                groups(k(br(i)+1:br(i+1)))=nGroups;
                dependencyGroups(nGroups,:)=dependencyGroups(g,:);
            end
        end
    end
    %[(1:nInstructions)',oldgroups,groups]
else
    %fprintf('WARNING: break of non-contiguous dependency groups has been bypassed, this could perhaps lead to cyclic dependencies\n')
end

%% Make sure there are no cyclic dependencies among dependency groups
nzgroups=find(groups>0);
in=1:nInstructions;
p=sparse(in(nzgroups),groups(nzgroups),ones(length(nzgroups),1),nInstructions,nGroups);
groupsGraph=p'*obj.dependencyGraph*p;
[i,j]=find(groupsGraph);
k=find(i==j);
i(k)=[];
j(k)=[];
groupsGraph=sparse(i,j,ones(size(i)),nGroups,nGroups);
%full(groupsGraph)
if 0
    % Use bioinformatics toolbox
    try
        computeOrder=graphtopoorder(groupsGraph');
    catch me
        error('dependencies graph appears to have cycles, this should not happen ;-(');
    end
else
    % Use Transitive reduction of a DAG by Frederik Gwinner
    computeOrder=toporder(groupsGraph');
    if isempty(computeOrder)
        error('dependencies graph appears to have cycles, this should not happen ;-(');
    end
end

if verboseLevel>2
    fprintf('Original ordering:\n');
    for i=1:length(computeOrder)
        fprintf('Group %d [%s]: %s\n',...
                i,index2str(dependencyGroups(i,:)),index2str(find(groups==i)));
    end

    fprintf('Before re-ordering:\n');
    for i=1:length(computeOrder)
        fprintf('Group %d [%s]: %s\n',computeOrder(i),index2str(dependencyGroups(computeOrder(i),:)),index2str(find(groups==computeOrder(i))));
    end
end

fprintf(' (%.2f sec), ordering groups',etime(clock,t0));

%% reorder groups to get the right computation order
dependencyGroups=dependencyGroups(computeOrder,:);
p=zeros(length(computeOrder),1);
p(computeOrder)=1:length(computeOrder);
k=find(groups>0);
groups(k)=p(groups(k));

if verboseLevel>2
    fprintf('After re-ordering:\n');
    for i=1:length(computeOrder)
        fprintf('Group %d [%s]: %s\n',i,index2str(dependencyGroups(i,:)),index2str(find(groups==i)));
    end
end

%% Mark sets/gets/copies with the appropriate dependency groups
j=1;
for i=1:length(obj.gets)
    obj.gets(i).parentGroups=find(dependencyGroups(:,j));
    j=j+1;
end
for i=1:length(obj.saves)
    obj.saves(i).parentGroups=find(dependencyGroups(:,j));
    j=j+1;
end
for i=1:length(obj.copies)
    obj.copies(i).parentGroups=find(dependencyGroups(:,j));
    j=j+1;
end
for i=1:length(obj.sets)
    obj.sets(i).childrenGroups=find(dependencyGroups(:,j));
    j=j+1;
end
for i=1:length(obj.copies)
    obj.copies(i).childrenGroups=find(dependencyGroups(:,j));
    j=j+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% determine input, output, and local variables for the different groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(' (%.2f sec), local variables',etime(clock,t0));

if true
    maxLocal=0;
    gT=obj.dependencyGraph';
    nonLocalInstructions=false(1,size(gT,1));
    for g=1:nGroups
        k=(groups==g);

        childrenOutsideGroup=obj.dependencyGraph(~k,k);
        parentsOutsideGroup=gT(~k,k);   % faster than using obj.dependencyGraph(k,~k);

        k=find(k);

        obj.groupInputInstructions{g}=union(k(any(parentsOutsideGroup,1)),...
                                            intersect(inputInstructions,k));
        obj.groupOutputInstructions{g}=union(k(any(childrenOutsideGroup,1)),...
                                            intersect(outputInstructions,k));
        obj.groupLocalInstructions{g}=setdiff(k,...
                                              union(obj.groupInputInstructions{g},...
                                                    obj.groupOutputInstructions{g}));
        maxLocal=max(maxLocal,length(obj.groupLocalInstructions{g}));
        nonLocalInstructions(obj.groupInputInstructions{g})=true;
        nonLocalInstructions(obj.groupOutputInstructions{g})=true;
        if verboseLevel>0
            fprintf(' group %3d [%6d..%6d] has %6d variables: input %6d, output %6d, local %6d (non local %6d/%6d, max local %6d)\n',...
                    g,min(k),max(k),length(k),...
                    length(obj.groupInputInstructions{g}),...
                    length(obj.groupOutputInstructions{g}),...
                    length(obj.groupLocalInstructions{g}),...
                    sum(nonLocalInstructions),length(nonLocalInstructions),maxLocal);
            %fprintf('        local instructions:\n');
            %disp(obj.groupLocalInstructions{g})
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute mapping between instructions and memory locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(' (%.2f sec), memory locations (',etime(clock,t0));

obj.memoryLocations=(1:nInstructions);

%obj.memoryLocations=randperm(nInstructions);

%obj.memoryLocations=(nInstructions:-1:1);

%obj.memoryLocations=(2:nInstructions+1);

if reuseLocalMemory
    fprintf(' reuse-across-groups');
    %% memory location reusing local variables
    nFixed=sum(nonLocalInstructions);
    % non-local instructions
    obj.memoryLocations(nonLocalInstructions)=1:nFixed;
    % local instructions (reuse from group to group)
    for g=1:nGroups
        obj.memoryLocations(obj.groupLocalInstructions{g})=...
            nFixed+(1:length(obj.groupLocalInstructions{g}));
    end
end

if findUnusedMemory
    fprintf(' unused-instructions');
    %% find unused instructions
    k=find((groups==0) & ~findInstructionsByType(obj.Itypes.I_set));
    nUnused=length(k);
    obj.memoryLocations(k)=NaN;

    if verboseLevel>0
        fprintf('\n  dependencyGroups: %d instructions not used\n',nUnused);
    end
end

nReuses=0;
nLocalVariables=0;
if reuseMemory
    fprintf(' reuse-within-group');
    %% look for variables internal to the group, and reuse memory when possible
    nGroups=size(dependencyGroups,1);
    for g=1:nGroups
        k=(groups==g);

        childrenOutsideGroup=obj.dependencyGraph(~k,k);
        k=find(k);

        % local variables used by group -- no children outside group
        localVariables=k(~any(childrenOutsideGroup,1));
        % and not outputs
        localVariables=setdiff(localVariables,outputInstructions);

        % last time each local variable is used
        if ~isempty(localVariables)
            usedBy=obj.dependencyGraph(:,localVariables);   % 1 if used
            [i,j]=find(usedBy);
            usedBy=sparse(i,j,i,size(usedBy,1),size(usedBy,2)); % variable # if used
            % fprintf('localVariables;usedBy\n');
            % disp(full([localVariables';usedBy]))
            lastUsedVariable=full(max(usedBy,[],1));
            % fprintf('localVariables;lastUsedVariable\n');
            % disp([localVariables';lastUsedVariable])
            lastUsedVariable(lastUsedVariable==0)=Inf;
            % memory used by each local variable
            localMemory=obj.memoryLocations(localVariables);
            lastUsedMemory=lastUsedVariable;
        end

        if verboseLevel>1
            fprintf(' group %3d [%6d..%6d] has %4d/%4d local variables',...
                    g,min(k),max(k),...
                    length(localVariables),length(k));
            if verboseLevel>2
                fprintf('\n');
            end
        end

        for j=1:length(localVariables)
            if verboseLevel>2
                fprintf('  instruction %d, last used by %d\n',...
                        localVariables(j),lastUsedVariable(j));
            end

            if verboseLevel>2
                fprintf('   localVariables;lastUsed;memoryLocations\n');
                disp([localVariables';
                      lastUsedVariable;
                      obj.memoryLocations(localVariables)]);
                fprintf('   localMemory;lastUsed;lastUsed<localVariables(%d)\n',j);
                disp([localMemory;
                      lastUsedMemory;
                      lastUsedMemory<localVariables(j)]);
            end
            cond=lastUsedMemory<localVariables(j);
            if any(cond)
                i=find(cond,1,'first');
                if verboseLevel>2
                    fprintf('     reusing memory location %d (instead of %d), for instruction %d\n',...
                            localMemory(i),...
                            obj.memoryLocations(localVariables(j)),...
                            localVariables(j));
                end
                % update memory
                obj.memoryLocations(localVariables(j))=localMemory(i);
                % update last time used
                lastUsedMemory(i)=lastUsedVariable(j);
                % remove unused memory from pool
                lastUsedMemory(j)=inf;  % prevent from being reused
                nReuses=nReuses+1;
            end
        end

        nLocalVariables=nLocalVariables+length(localVariables);
        if verboseLevel>1
            if verboseLevel>2
                fprintf('   ');
            else
                fprintf(', ');
            end
            fprintf('  re-used %4d memory locations (%4d local variables)\n',...
                    nReuses,nLocalVariables);
        end
    end
end

fprintf(' %.2f sec) ',etime(clock,t0));

if verboseLevel>0
    fprintf('  dependencyGroups: re-used %d memory locations (%4d local variables)\n',...
            nReuses,nLocalVariables);
end

if compressMemory
    %% remove memory holes
    if verboseLevel>2
        obj.memoryLocations
    end
    if verboseLevel>0
        fprintf('  dependencyGroups: compressed memory from %d',...
                max(obj.memoryLocations));
    end
    [~,~,obj.memoryLocations]=unique(obj.memoryLocations,'stable');
    obj.memoryLocations=obj.memoryLocations(:)';
    if verboseLevel>2
        obj.memoryLocations
    end
    if verboseLevel>0
        fprintf(' to %d\n',max(obj.memoryLocations));
    end
end

%% Save groups in object
obj.dependencyGroups=dependencyGroups;
%set(obj.instructions,'group',1:nInstructions,groups);
obj.instructionsGroup=groups;

obj.outputInstructions=outputInstructions;
obj.inputInstructions=inputInstructions;
obj.nonLocalInstructions=nonLocalInstructions;

if verboseLevel>1
    fprintf('Dependency groups:\n');
    fprintf('%36s','group');
    fprintf('%2d',1:size(obj.dependencyGroups,1));
    fprintf('\n');
    fprintf('%s\n',repmat('-',1,35));
    for i=1:size(obj.dependencyGroups,2)
        fprintf('%-36s',obj.dependencyGroupColName{i});
        fprintf('%2d',full(obj.dependencyGroups(:,i)));
        fprintf('\n');
    end
    %disp(obj,1);
end

%% collect statistics
obj.statistics.nGroups=nGroups;
obj.statistics.nSets=length(obj.sets);
obj.statistics.nGets=length(obj.gets);
obj.statistics.nCopies=length(obj.copies);
obj.statistics.nInstructions=nInstructions;
obj.statistics.sizeScratchbook=max(obj.memoryLocations);

end