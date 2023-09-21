function [KLS, P, Q] = nk_KLS2(S,T)
% =========================================================================
% function [KLS, P, Q] = nk_KLS(S,T)
% =========================================================================
% This functions computes the symmetric KL divergence between two
% one-dimensional distributions following the descriptions of Wang et al.
% Brain and Behaviour, 2016, 6(4): e00448
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 03/2021

% Compute kernel-density estimates
[~,Sd] = kde(S); [~,Pd] = kde(T);
% Compute PDFs over the density estimates
P = normpdf(Sd,mean(Sd),std(Sd))+eps; Q = normpdf(Pd,mean(Pd),std(Pd))+eps;
% Compute symmetric KL divergence
KLS = exp(-sum(P.*log(P./Q) + Q.*log(Q./P)));

