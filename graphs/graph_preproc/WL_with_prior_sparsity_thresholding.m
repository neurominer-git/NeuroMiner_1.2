function K = WL_with_prior_sparsity_thresholding(Y1, Y2, sp, param1) % in case of test kernel, Y1 should be Ytest
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
        spArray = apply_sparsity_thres2(Y(i,:), sp);
        % transform to 2D matrix
        graph = array_to_graph(spArray); 
        
        gList = [gList, graph]; 
        
    end
    lab = 0; 
    KM = feval("WL", gList, param1, lab);
    K = KM{size(KM,2)}; % which one to choose?
    if test % we want test kernel of dimension n_test x n_train
        K = K(1:size(Y1,1),1:size(Y2,1));
    else
        K = [(1:size(K,1))', normalizekm(K)];
    end     
end