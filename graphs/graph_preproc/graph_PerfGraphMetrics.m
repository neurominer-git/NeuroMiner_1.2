function [sY,IN] = graph_PerfGraphMetrics(Y,IN)

% =========================== WRAPPER FUNCTION ============================
if ~exist('IN','var'), IN = []; end
if iscell(Y)
    sY = cell(1,numel(Y));
    for i=1:numel(Y), [sY{i}, IN] =  PerfGraphMetrics(Y{i}, IN); end
else
    [ sY, IN ] = PerfGraphMetrics(Y, IN);
end
end

%-------------------------------------------------------------------------
function [Y, IN] = PerfGraphMetrics(Y, IN)
netwmetrics{1}.id           = 'degree';
netwmetrics{1}.listidx      =  1;
netwmetrics{1}.func         = 'degree = degrees_und(A);';
netwmetrics{2}.id           = 'strength';
netwmetrics{2}.listidx      =  2;
netwmetrics{2}.func         = 'strength = strengths_und(A);';
netwmetrics{3}.id           = 'betweenness';
netwmetrics{3}.listidx      =  3;
netwmetrics{3}.func         = 'betweenness = betweenness_wei(A);';
netwmetrics{4}.id           = 'clustering_coef';
netwmetrics{4}.listidx      =  4;
netwmetrics{4}.func         = 'clustering_coef = clustering_coef_wu(A);';
netwmetrics{5}.id           = 'clustcoef';
netwmetrics{5}.listidx      =  5;
netwmetrics{5}.func         = 'A(A >0) = 1; clustcoef = clustering_coef_bu(A);';
netwmetrics{6}.id           = 'diameter';
netwmetrics{6}.listidx      =  6;
netwmetrics{6}.func         = '[lambda,efficiency,ecc,radius,diameter] = charpath(A);';
% netwmetrics{7}.id           = 'small-worldness (global)';
% netwmetrics{7}.listidx      =  7;
% netwmetrics{7}.func         = @degrees_und;
netwmetrics{7}.id           = 'transitivity';
netwmetrics{7}.listidx      =  7;
netwmetrics{7}.func         = 'transitivity = transitivity_wu(A);';
netwmetrics{8}.id           = 'eigenvec';
netwmetrics{8}.listidx      =  8;
netwmetrics{8}.func         = 'eigenvec = eigenvector_centrality_und(A);';
netwmetrics{9}.id           = 'efficiency_loc';
netwmetrics{9}.listidx      =  9;
netwmetrics{9}.func         = 'efficiency_loc = efficiency_wei(A,2);';
netwmetrics{10}.id          = 'efficiency';
netwmetrics{10}.listidx     =  10;
netwmetrics{10}.func        = 'efficiency = efficiency_wei(A);';
netwmetrics{11}.id          = 'closeness';
netwmetrics{11}.listidx     =  11;
netwmetrics{11}.func        = 'closeness = centrality(G, ''closeness'', ''Cost'', G.Edges.Weight);';
netwmetrics{12}.id          = 'radius';
netwmetrics{12}.listidx     =  12;
netwmetrics{12}.func        = '[lambda,eff,ecc,radius,diameter] = charpath(A)';
netwmetrics{13}.id          = 'pagerank';
netwmetrics{13}.listidx     =  13;
netwmetrics{13}.func        = 'pagerank = pagerank_centrality(A)';
netwmetrics{14}.id          = 'distance';
netwmetrics{14}.listidx     =  14;
netwmetrics{14}.func        = 'distance = distance_wei(A);';

local_metrics               = [1, 2, 3, 4, 8, 9, 11, 13, 14];
global_metrics              = [5, 6, 7, 10, 12];

% allocate space
% count local network metrics
sel_metrics                 = IN.metricslist;
locm = 0; 
globm = 0; 
funclist = {};
for i = 1:length(sel_metrics)
    thismetric = sel_metrics(i);
    thismetric = thismetric{1};
    if ismember(thismetric.listidx, local_metrics)
        locm = locm+1;
    else
        globm = globm+1;
    end

    funclist{i} = netwmetrics{thismetric.listidx}.func;

    if i == 1
        metvec_str = sprintf('metrics_vector = %s', netwmetrics{thismetric.listidx}.id);
    else
        metvec_str = sprintf('%s, %s', metvec_str, netwmetrics{thismetric.listidx}.id);
    end
    %varnamelist{i} = netwmetrics{sel_metrics(i).listidx}.id;
end

n_nodes = (1+sqrt(1+(8*size(Y,2))))/2; 
metricsY = zeros(size(Y,1), (n_nodes*locm + globm));%zeros(size(Y,1),(n_nodes*6+3));
for i = 1:size(Y,1)
    G_struct = array_to_graph(Y(i,:));
    %G = graph(G_struct.am);
    A = G_struct.am;

    for j = 1:length(sel_metrics)
        eval(funclist{j});
        eval(sprintf('%s;', metvec_str));
        %         %local
        %         degree = degrees_und(A);
        %         strength = strengths_und(A);
        %         betweenness = betweenness_wei(A);
        %         pagerank = pagerank_centrality(A);
        %         clustering_coef = clustering_coef_wu(A);
        %         transitivity = transitivity_wu(A);
        %         efficiency_loc = efficiency_wei(A,2);
        %         distance = distance_wei(A);
        %         [lambda,efficiency,ecc,radius,diameter] = charpath(A);
        %         efficiency = efficiency_wei(A);
        %         eigenvec = eigenvector_centrality_und(A);
        %         %closeness = centrality(G, 'closeness', 'Cost', G.Edges.Weight);
        %         %strength = centrality(G, 'degree', 'Importance',G.Edges.Weight);
        %         %degree = centrality(G, 'degree');
        %         %betweenness = centrality(G, 'betweenness', 'Cost', G.Edges.Weight);
        %         %eigenvector = centrality(G, 'eigenvector', 'Importance', G.Edges.Weight);
        %
        % %         A = G_struct.am;
        % %         A(A >0) = 1;
        % %         clustcoef = clustering_coef_bu(A);
        %
        %         % smallworldness: https://github.com/mdhumphries/SmallWorldNess/blob/master/Computing_And_Testing_S_On_Data_Network.m
        %         % get its basic properties
        % %         n = size(A,1);  % number of nodes
        % %         k = sum(A);  % degree distribution of undirected network
        % %         m = sum(k)/2;
        % %         K = mean(k); % mean degree of network
        % %         FLAG_Cws = 1;
        % %         [expectedC,expectedL] = ER_Expected_L_C(K,n);  % L_rand and C_rand
        % %         [S_ws,C_ws,L] = small_world_ness(A,expectedL,expectedC,FLAG_Cws);  % Using WS clustering coefficient
        %
        % %         metrics_vector = [closeness', ...
        % %             strength', ...
        % % %             degree', ...
        % % %             betweenness', ...
        % % %             eigenvector', ...
        % %             clustcoef', ...
        %             S_ws, C_ws, L];
        %         metrics_vector = [degree, ...
        %             strength, ...
        %             betweenness', ...
        %             clustering_coef', ...
        %             eigenvec', ...
        %             efficiency_loc', ...
        %             ecc', ...
        %             efficiency, ...
        %             diameter, ...
        %             transitivity];


        metricsY(i,:) = metrics_vector;
    end
    Y = metricsY;
end
end


