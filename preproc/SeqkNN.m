function [imputed_data]=SeqkNN(data, K, traindata)
% SeqKNN: Sequential KNN imputation method
% This function estimates missing values sequentially from the gene that has
% least missing rate in microarray data, using weighted mean of k nearest neighbors.
%
% <Usage>
% imputed_data=SeqKNN(data, k);
%
% <Arguments>
% data: matrix or dataframe, 1 row corresponds to 1 gene, 1 column to 1
% sample,colnames and rownames can be used
% k: number of nearest neighbors
%
% <Details>
% SeqKNN separates the dataset into incomplete and complete sets that have
% or have not missing values, respectively. 
% The genes in incomplete set are imputed by the order of missing rate. Missing value
% is filled by the weighted mean value of corresponding column of the nearest neighbor genes in
% complete set. Once all missing values in a gene are imputed, the imputed gene is moved into the
% complete set and used for the imputation of the rest of genes in incomplete set. In this process,
% all missing values in one gene can be imputed simultaneously from the selected neighbor genes
% in complete set. This reduces execution time from previously developed KNN method that selects
% nearest neighbors for each imputation.
%
% <References>
% Ki-Yeol Kim, Byoung-Jin Kim, Gwan-Su Yi (2004.Oct.26) "Reuse of imputed data in microarray
% analysis increases imputation efficiency", BMC Bioinformatics 5:160.
% -------------------------------------------------------------------------
%%
imputed_data=zeros(size(data));
complete=[];
incomplete=[];
missing=[];
com_ind=[];
incom_ind=[];
[rows,~]=size(data);

for i=1:rows
    if ~isnan(sum(data(i,:)))
        complete=[complete; data(i,:)];
        com_ind=[com_ind i];
    else
        incomplete=[incomplete; data(i,:)];
        incom_ind=[incom_ind i];
        missing=[missing sum(isnan(data(i,:)))];
    end
end   
imputed_data(com_ind,:)=complete;
[~,missing_order]=sort(missing);

if exist("traindata","var") && ~isequaln(data,traindata)
    completedata = traindata;
else
    completedata = complete;
end

%% part2
[irows,icols]=size(incomplete);
for j=1:irows
    fprintf('.');
    dist=[];
    cgen=size(completedata,1);
    for i=1:cgen
        dist(i)=nm_nansum((incomplete(missing_order(j),:)-completedata(i,:)).^2);
    end
    [dist,pos]=sort(dist);
    pos=pos(1:K);
    dist=dist(1:K);
    weights=(1./dist)./sum(1./dist);
    for g=1:icols
        if (isnan(incomplete(missing_order(j),g)))
            incomplete(missing_order(j),g)=weights*(completedata(pos,g));
        end
    end
    completedata = [completedata; incomplete(missing_order(j),:)];
end
imputed_data(incom_ind,:)=incomplete;
end