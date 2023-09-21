function C = clusttriang(A)

% CLUSTTRIANG computes clustering coefficient based on number of triangles
%   C = CLUSTTRIANG(A) compute the clustering coefficient C of the
%   adjacency matrix A, based on the ratio (number of triangles)/(number of triples).
%   
%   Notes:
%   (1)This is a subtly different definition than that of Watts & Strogatz (1998). 
%   In the present form, it does give the fraction of neighbouring nodes
%   that are also neighbours of each other. (Also note: for directed
%   graphs, this means *any* direction of connection, as along as all three
%   nodes are connected).
%
%   (2) This function uses the method of Keeling (1999) to compute C. 
%
%   Mark Humphries 22/8/2006

A2 = A * A;
A3 = A2 * A;
sumA2 = sum(sum(A2));

C = trace(A3) / (sumA2 - trace(A2));