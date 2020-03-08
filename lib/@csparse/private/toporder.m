function node_order = toporder(adjmat)
% node_order = toporder(adjmat)
%
%   Computes a topological ordering of the nodes for a directed
%   acyclic graph defined by its adjacency matrix.
%
% Adapted from:
%    Transitive reduction of a DAG by Frederik Gwinner
%    Copyright (c) 2011, Frederik Gwinner
%    All rights reserved.
%    Code covered by the BSD License

    % Check whether adjacency matrix is square
    V=size(adjmat,1);
    
    % % get topological ordering of the nodes
    node_indices = 1:V;
    node_order=[];
    start_nodes = node_indices(sum(adjmat,1)==0);%nodes having no incoming edges
    sort_adj=adjmat; % copy of the initial matrix
    while ~isempty(start_nodes)
        n = start_nodes(1);
        start_nodes(1)=[];
        node_order=[node_order,n];
        successors = node_indices(logical(adjmat(n,:)));
        for mi=1:numel(successors)
            m=successors(mi);
            sort_adj(n,m)=0;
            if sum(sort_adj(:,m))==0
                start_nodes=[start_nodes,m];
            end
        end
    end

    if sum(sum(sort_adj))~=0
        % if there's still an edge left, the graph must have
        % contained a cycle
        node_order=[];
    end

end