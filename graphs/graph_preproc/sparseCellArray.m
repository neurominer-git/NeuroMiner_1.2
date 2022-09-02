function gCell = sparseCellArray(Y1, Y2, sp) % in case of test kernel, Y1 should be Ytest
    if ~isequal(Y1,Y2)
        test = 1;
        Y = vertcat(Y1,Y2);
    else
        test = 0;
        Y = Y1;
    end
    gCell= {};
    for i = 1:size(Y,1) 
        % apply sparsity threshold
        spArray = apply_sparsity_thres(Y(i,:), sp);
        % transform to 2D matrix
        graph = array_to_graph2(spArray); 
        gCell(1,i) = {graph};
        %gList = [gList, graph]; 
        
    end
end