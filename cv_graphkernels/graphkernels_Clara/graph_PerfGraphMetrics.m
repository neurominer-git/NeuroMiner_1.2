function [sY,IN] = graph_PerfGraphMetrics(Y,IN)

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfGraphMetrics(Y{i}, IN); end
    else
        [ sY, IN ] = PerfGraphMetrics(Y, IN );
    end
end
function [Y, IN] = PerfGraphMetrics(Y, IN)
    n_nodes = (1+sqrt(1+(8*size(Y,2))))/2;
    metricsY = zeros(size(Y,1), (n_nodes*7+3));%zeros(size(Y,1),(n_nodes*6+3));
    for i = 1:size(Y,1)
        G_struct = array_to_graph(Y(i,:));
        %G = graph(G_struct.am);
        A = G_struct.am;
        %local
        degree = degrees_und(A);
        strength = strengths_und(A);
        betweenness = betweenness_wei(A);
        %pagerank = pagerank_centrality(A);
        clustering_coef = clustering_coef_wu(A);
        transitivity = transitivity_wu(A);
        efficiency_loc = efficiency_wei(A,2);
        %distance = distance_wei(A);
        [lambda,efficiency,ecc,radius,diameter] = charpath(A);
        efficiency = efficiency_wei(A);
        eigenvec = eigenvector_centrality_und(A);
        %closeness = centrality(G, 'closeness', 'Cost', G.Edges.Weight);
        %strength = centrality(G, 'degree', 'Importance',G.Edges.Weight); 
        %degree = centrality(G, 'degree');
        %betweenness = centrality(G, 'betweenness', 'Cost', G.Edges.Weight);
        %eigenvector = centrality(G, 'eigenvector', 'Importance', G.Edges.Weight);
        
%         A = G_struct.am;
%         A(A >0) = 1; 
%         clustcoef = clustering_coef_bu(A);
        
        % smallworldness: https://github.com/mdhumphries/SmallWorldNess/blob/master/Computing_And_Testing_S_On_Data_Network.m
        % get its basic properties
%         n = size(A,1);  % number of nodes
%         k = sum(A);  % degree distribution of undirected network
%         m = sum(k)/2;
%         K = mean(k); % mean degree of network
%         FLAG_Cws = 1;
%         [expectedC,expectedL] = ER_Expected_L_C(K,n);  % L_rand and C_rand
%         [S_ws,C_ws,L] = small_world_ness(A,expectedL,expectedC,FLAG_Cws);  % Using WS clustering coefficient
        
%         metrics_vector = [closeness', ...
%             strength', ...
% %             degree', ...
% %             betweenness', ...
% %             eigenvector', ...
% %             clustcoef', ...
%             S_ws, C_ws, L];
        metrics_vector = [degree, ...
            strength, ...
            betweenness', ...
            clustering_coef', ...
            eigenvec', ...
            efficiency_loc', ...
            ecc', ...
            efficiency, ...
            diameter, ...
            transitivity];
            
            
        metricsY(i,:) = metrics_vector;
    end
Y = metricsY;
end


