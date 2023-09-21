function A = random_graph(N,I,P,type,varargin)

% RANDOM_GRAPH create a random graph
%
%   A = RANDOM_GRAPH(N,I,P,T) where N is the number of vertices in
%   the graph, I is the proportion of inhibitory vertices (useful for asymmetry or sign problems),
%   P is the connection probability, and T is the type of edges that 
%   constitute the graph ('directed' or 'undirected'). The graph created is entirely random,
%   with the probability of connection tested for each pair of vertices
%   (Erdos-Renyi random graph). Only single edges and no self-loops are
%   allowed - hence, P is adjusted to account for the slight reduction in
%   the number of maximum edges.
%
%   If P > 1 then it is assumed that P is the number of unique edges to be made. That is, the number of edges is
%   specified and connect randomly selected pairs - this allows the creation of 'equivalent' random
%   graphs where the average probability P(k) (connections per node) is unknown or meaningless (such as when an equivalent random
%   graph is required for a structure with a double Gaussian distribution of P(k)). Note that
%   in this case P < N(N-1)/2 for undirected graphs, and P < N(N-1) for
%   directed graphs. 
%
%   A = RANDOM_GRAPH(...,S) initialises the random number generator with seed S, allowing for repeated
%   creation of same graph (set = [] to omit) 
%
%   Returns an adjacency matrix A which defines every individual connection (ROW is connections
%   from the neuron, COLUMN connections to the neuron); for undirected graphs, the connections are
%   symmetrical about the diagonal. Entries are {0,1,-1}. 
%
%   Mark Humphries 29/11/2006

if ~strcmp(lower(type),'directed') & ~strcmp(lower(type),'undirected')
    error('Incorrect graph type specified')
end

% set seed for random number generator
if nargin >= 5
    S = varargin{1};
else
    S = [];    
end

    
if P <= 1    
    % set random number generator
    if S
        rand('state',S);
    end

	P = (N^2 * P) / (N^2 - N);     %% re-calculate connection probability to account for removal of self-connections below
                                               %% note: this approximation fails as P -> 1 for small N    
	%%%%% adjacency matrix
    if strcmp(type,'undirected')
        A = triu(rand(N));
    else
        A = rand(N);
    end
  
    A(A > 0 & A <= P) = 1;                 
    A(A < 1) = 0;
    A(eye(N)==1) = 0;
    
else    
    if strcmp(type,'undirected')
        if P >= N*(N-1)/2
             error('Too many edges for the specified number of nodes in an undirected graph');
        end
        link_pairs = rand_pair(1:N,1:N,P,S,'bi','diff');                       % only unique node pairs allowed for undirected  
    else
        if P >= N*(N-1)
             error('Too many edges for the specified number of nodes in a directed graph');
        end
        link_pairs = rand_pair(1:N,1:N,P,S,[],'diff');                         % any pairs allowed
    end

    link_ind = sub2ind([N,N],link_pairs(:,1),link_pairs(:,2));                 % convert to linear index for use in matrix
    
    A = zeros(N,N);
    A(link_ind) = 1;
end
                                                                
if strcmp(type,'undirected')
    % duplicate rows in columns (so that all out connections are also in connections)
    A = A + A';
    %A(A>1) = 1;     % if already symmetric then reset to 1!
end

%%% determine which neurons are inhibitory
num_I = round(N*I); % number of inhibitory neurons
I_indices = randperm(N);
I_indices = I_indices(1:num_I);

%%%%% weight matrix
A(I_indices,:) = A(I_indices,:) .* -1; 
if strcmp(type,'undirected')
    A(:,I_indices) = A(:,I_indices) .* -1;
end

