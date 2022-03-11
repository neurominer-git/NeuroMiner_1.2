function [Lrand,Crand] = NullModel_L_C(n,m,Nrepeats,FLAG)

% NULLMODEL_L_C Monte Carlo estimates of path-length and clustering of an ER random graph
% [C,L] = NULLMODEL_L_C(N,M,R,FLAG) for N nodes and M edges,
% creates R realisations of an Erdos-Renyi random graph (N,M), and
% computes the shortest path-length L and clustering coefficient C 
% for each one
%
% FLAG is a number indicating which clustering coefficient to compute:
%   1 - Cws 
%   2 - transitivity C (no. of triangles)

% Mark Humphries 3/2/2017


for iE = 1:Nrepeats
    ER = random_graph(n,0,m,'undirected');  % make E-R random graph
    [~,D] = reachdist(ER);  % returns Distance matrix of all pairwise distances
    Lrand(iE) = mean(D(:));  % mean shortest path-length: including self-loops
    
%     if isinf(Lrand(iE))
%         keyboard
%     end
    
    % calculate required form of C
    switch FLAG
        case 1
            c = clustering_coef_bu(ER);  % vector of each node's C_ws
            Crand(iE) = mean(c);  % mean C
        case 2
            Crand(iE) = clusttriang(ER);
    end
end