function graph = array_to_graph(a)

NODES = (1+sqrt(1+(8*length(a))))/2;
adj_m = zeros(NODES);

%a(a>0) = 1;
adj_m(triu(true(NODES),1)) = a;
adj_m = adj_m + transpose(adj_m);
graph.am = adj_m; %reshape(a,[sqrt(length(a)),sqrt(length(a))]);
graph.al = cellfun(@(x) find(x),num2cell(graph.am,2),'un',0);
nodes = 1:size(adj_m,1);
graph.nl.values = nodes';
end


