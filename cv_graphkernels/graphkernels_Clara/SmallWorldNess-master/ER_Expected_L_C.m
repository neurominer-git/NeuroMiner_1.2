function [expectedC,expectedL] = ER_Expected_L_C(k,n)
    
% ER_EXPECTED_L_C the expected path-length and clustering of an ER random graph
% [C,L] = ER_EXPECTED_L_C(K,N) for a network of N nodes and mena degree K, 
% computes the expected shortest path-length L and expected clustering 
% coefficient C if that network was an Erdos-Renyi random graph
%
% Mark Humphries 3/2/2017

expectedC = k / n;
z1 = k;
z2 = k^2;
expectedL = (log((n-1)*(z2 - z1) + z1^2) - log(z1^2)) / log(z2/z1);
