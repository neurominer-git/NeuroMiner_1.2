
function graph = array_to_graph2(a)
% binarize adjacency matrix
a(a>0) = 1;

am = reshape(a,[sqrt(length(a)),sqrt(length(a))]);
graph.edges = adj2edge(am);
%graph.edges = cellfun(@(x) find(x),num2cell(am,2),'un',0);
nodes = 1:sqrt(length(a));
graph.nodelabels = nodes';
end


function el=adj2edge(adj)
    n=length(adj); % number of nodes
    edges=find(triu(adj>0));  % indices of all edges
    el=[];
    for e=1:length(edges)
        [i,j]=ind2sub([n,n],edges(e)); % node indices of edge e  
        el=[el; i j adj(i,j)];
    end
    el = el(:,1:2);
end