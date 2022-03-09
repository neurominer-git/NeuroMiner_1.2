function [alpha,P,Sb] = SsampleER(n,k,m,I,ns,S,FLAG)

% SSAMPLEER generate small-world-ness values for E-R random graph
% [A,P,V] = SSAMPLEER(N,K,M,I,SAMP,S,FLAG) computes SAMP realisations of an
%   E-R random graph, of size N, mean degree K, and number of unique edges M. 
%   For each, the small-world-ness is computed. 
%   Returned:
%       A: the I% interval of the null model Small-world-ness distribution.  I should be given as proportion
%           i.e. I = 0.95 is the 95% interval. This measures the extremes in variation of S under the null model 
%       P: the P-value for the one-sided test for the hypothesis that the data small-world-ness S is 
%           greater than expected from the null-model small-world-ness distribution.
%           (Technically, that S is incorrectly rejected from the null-model distirbution) 
%       V: the vector of all null model small-world-ness values computed
%
%   The value of FLAG determines which clustering coefficient is computed:
%       1 - raw S: clustering coefficient of Watts & Strogatz (1998)
%       2 - raw S: transitivity clustering coefficient (no. of triangles) e.g.
%           Barat & Weigt(2000); Newman et al (2001);   
%
%   Notes:
%   (1) The S-score computed here is based on comparison with the E-R random graph, following the
%       Watts & Strogatz (1998) definition.
%   (2) The expected C and L values for an E-R random graph are computed, implying:
%       (i) Small n would lead to errors, as path length is underestimated 
%       (ii) This procedure is currently only defined for undirected graphs,
%
%   Mark Humphries 3/2/2017

[expectedC,expectedL] = ER_Expected_L_C(k,n);

for i = 1:ns
    A = random_graph(n,0,m,'undirected');
    [SER(i),~,~] = small_world_ness(A,expectedL,expectedC,FLAG);
end

% remove any SER = 0, due to INF path-length
Sb = SER(SER > 0);
Nsamps = numel(Sb);

% confidence interval on the null model value of S
sig = (1 - I);


% keyboard

alpha = prctile(Sb,100*[sig/2 1-sig/2]);

% P-value for S
P = sum(Sb > S)+1/Nsamps;




