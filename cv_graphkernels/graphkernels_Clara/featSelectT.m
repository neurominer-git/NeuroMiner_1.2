function selectedY = featSelectT(X1, X2, y, train) % y = labels

    if train
        
        g1_idx = find(y==-1);
        g2_idx = find(y==1);
        selectedY = X1; 
        for i = 1:size(X1,2)
            g1 = X1(g1_idx,i);
            g2 = X1(g2_idx,i);
            [h,p] = ttest2(g1,g2);
            if h == 0 || isnan(h)
                selectedY(:,i) = 0;
            else 
                i
                h
            end
                
        end
    else
        selectedY = X2;
        for i = 1:size(X2,2)
            if sum(X1(:,i) == 0
                selectedY(:,i) = 0;
            end
        end
    end
end
