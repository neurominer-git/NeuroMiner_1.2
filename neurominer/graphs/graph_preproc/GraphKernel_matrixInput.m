function K = GraphKernel_matrixInput(Y1, Y2, kernelf, param1) % in case of test kernel, Y1 should be Ytest
    if ~isequal(Y1,Y2)
        test = 1;
        Y = vertcat(Y1,Y2);
    else
        test = 0;
        Y = Y1;
    end
    gList = [];
    for i = 1:size(Y,1) 
        % apply sparsity threshold
        % spArray = apply_sparsity_thres(Y(i,:), sp);
        % transform to 2D matrix
        graph = array_to_graph(Y(i,:)); 
        gList = [gList, graph]; 
    end
    lab = 1; 
    KM = feval(kernelf, gList, param1, lab);
    K = KM{size(KM,2)}; % which one to choose?
    if test % we want test kernel of dimension n_test x n_train
        K = K(1:size(Y1,1),size(Y1,1)+1:size(Y1,1) + size(Y2,1));
    else
        K = [(1:size(K,1))', normalizekm(K)];
    end     
end