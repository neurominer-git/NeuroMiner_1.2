function [S,lambda,gamma] = SanalyticER(n,k,C,L,FLAG)

% SANALYTICER computes small-world-ness for analytic E-R graph values
%   [S,lambda,gamma] = SANALYTICER(N,K,C,L,FLAG) computes small-world-ness S, for networks of
%   size N, mean degree K, clustering coefficient C, and path length L,
%   where each is an array of values. The S score is based on the 
%   corresponding C and L values for the equivalent
%   Erdos-Renyi random graph (computed analytically). The form of S
%   computed depends on the value of FLAG:
%       1 - raw form (just S = (lambda / gamma))
%       2 - magnitude form (S = log10(lambda / gamma))
%
%   Will also return the path length ratio (lambda), and the clustering
%   coefficient ratio (gamma).
%
%   Note#1: Whether the resulting S value is for S(tri) or S(ws) depends
%   entirely on which C value is supplied by the user - the E-R graph
%   value is equivalent for either C(tri) or C(ws)
%
%   Note#2: for directed E-R random graphs, k should be mean in-degree; note that this
%   analytic form more susceptible to error as it relies on a giant component 
%   covering all of the graph, with very small in- and out-only populations.	
%
%   Mark Humphries 30/2/2007

C_ER = k ./ n;
L_ER = log(n) ./ log(k);

lambda = L ./ L_ER;
gamma = C ./ C_ER;

S = gamma ./ lambda;
if FLAG == 2
    S = log10(S);
end

