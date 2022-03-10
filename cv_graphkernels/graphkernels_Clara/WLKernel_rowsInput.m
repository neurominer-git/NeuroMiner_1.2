function K = WLKernel_rowsInput(Y, height)
    %K = zeros(size(Y,1));
    %size(K)
    gList = []
    for i = 1:size(Y,1)
        i
        %for j = 1:size(Y,1)
            
        graph = array_to_graph(Y(i,:)); 
            %gB = array_to_graph(Y(j,:)); 
        gList = [gList, graph]; 
    end
    %height = 3;
    lab = 0; 
    KM = WL(gList, height, lab);
    K = KM{size(KM,1)}; % which one to choose?
        
end
