function K = WL_with_prior_sparsityThres_Ttest(y, Y1, Y2, sp, param1) % in case of test kernel, Y1 should be Ytest
    
    if ~isequal(Y1,Y2)
        test = 1;
        spY1 = Y1;
        for i = 1:size(Y1,1)
            % apply sparsity threshold
            spArray = apply_sparsity_thres(Y1(i,:), sp);
            spY1(i,:) = spArray;
        end
        spY2 = Y2; 
        for i = 1:size(Y2,1)
            spArray = apply_sparsity_thres(Y2(i,:),sp);
            spY2(i,:) = spArray;
        end
        % feature selection with t-test
        selectedY1 = featSelectT(spY1, spY1, y, 1); % training
        selectedY1 = featSelectT(spY1, spY2, y, 0); % testing
        Y = vertcat(selectedY1,selectedY2);
    else
        test = 0;
        spY1 = Y1;
        for i = 1:size(Y1,1)
            % apply sparsity threshold
            spArray = apply_sparsity_thres(Y1(i,:), sp);
            spY1(i,:) = spArray;
        end
        selectedY1 = featSelectT(spY1, spY1, y, 1); % training
        Y = selectedY1;
    end
    %gCell = {};
    gList = [];
    for i = 1:size(Y,1)
        graph = array_to_graph(Y(i,:)); 
        %gCell(1,i) = {graph};
        gList = [gList, graph];
    end
    lab = 0; 
    KM = feval("WL", gList, param1, lab);
    KM
    K = KM{size(KM,2)}; % which one to choose?
    
    if test % we want test kernel of dimension n_test x n_train
        K = K(1:size(Y1,1),1:size(Y2,1));
    else
        K = [(1:size(K,1))', normalizekm(K)];
    
    end     
end