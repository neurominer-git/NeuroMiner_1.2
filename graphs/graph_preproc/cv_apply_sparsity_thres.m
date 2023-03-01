function S = cv_apply_sparsity_thres(A, thres)
% if mod(sqrt(size(A,2),1)) == 0 % the full connectivity matrix is entered
%     nEdges = (size(A,2)*thres)-size(A,2); % number of edges minus diagonal
% else % only upper triangle (without diagonal was entered)
nEdges = size(A,2)*thres; % number of edges
% end
nEdges = ceil(nEdges);
%sorteA = sort(A, 'descend');

S = A;
for i = 1:size(A,1)
    row = A(i,:);
    [B,I] = maxk(row,nEdges);
    sparseRow = row; 
    sparseRow(setdiff(1:length(row),I)) = 0;
    S(i,:) = sparseRow;
end
end
