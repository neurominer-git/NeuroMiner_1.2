%% script to compute small-world-ness and do statistical testing on example data network
% Mark Humphries 3/2/2017

clear all; close all;

% analysis parameters
Num_ER_repeats = 100;  % to estimate C and L numerically for E-R random graph
Num_S_repeats = 1000; % to get P-value for S; min P = 0.001 for 1000 samples
I = 0.95;

FLAG_Cws = 1;
FLAG_Ctransitive = 2;

%% load the adjacency matrix for the Lusseau bottle-nose dolphin social
% network
load dolphins  % loads struct of data in "Problem"; adjacency matrix is Problem.A

A = full(Problem.A); % convert into full from sparse format

% get its basic properties
n = size(A,1);  % number of nodes
k = sum(A);  % degree distribution of undirected network
m = sum(k)/2;
K = mean(k); % mean degree of network

%% computing small-world-ness using the analytical approximations for the E-R graph

[expectedC,expectedL] = ER_Expected_L_C(K,n);  % L_rand and C_rand

[S_ws,C_ws,L] = small_world_ness(A,expectedL,expectedC,FLAG_Cws);  % Using WS clustering coefficient
[S_trans,C_trans,L] = small_world_ness(A,expectedL,expectedC,FLAG_Ctransitive);  %  Using transitive clustering coefficient



%% computing small-world-ness by estimating L_rand and C_rand from an ensemble of random graphs
% check when using small networks...

[Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
[Lrand,CrandTrans] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Ctransitive);

% Note: if using a different random graph null model, e.g. the
% configuration model, then use this form

% compute small-world-ness using mean value over Monte-Carlo realisations

% NB: some path lengths in L will be INF if the ER network was not fully
% connected: we disregard these here as the dolphin network is fully
% connected.
Lrand_mean = mean(Lrand(Lrand < inf));

[S_ws_MC,~,~] = small_world_ness(A,Lrand_mean,mean(CrandWS),FLAG_Cws);  % Using WS clustering coefficient
[S_trans_MC,~,~] = small_world_ness(A,Lrand_mean,mean(CrandTrans),FLAG_Ctransitive);  %  Using transitive clustering coefficient


%% do Monte-Carlo test of disitribution of S in ER random graph

[I,P,Sb] = SsampleER(n,K,m,I,Num_S_repeats,S_ws,FLAG_Cws);  % for Sws here

% check how many samples ended up being used to calculate P-value
Nsamps = numel(Sb);
Pmax = 1 / Nsamps;