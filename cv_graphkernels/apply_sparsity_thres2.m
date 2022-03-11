function S = apply_sparsity_thres2(A, thres)
nEdges = size(A,2)*thres;
nEdges = ceil(nEdges);
%sorteA = sort(A, 'descend');

[B,I] = maxk(A,nEdges);
S = A;
S(setdiff(1:length(A),I)) = 0;
end