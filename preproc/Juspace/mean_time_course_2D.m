function [D,Reg_all] = mean_time_course_2D(data_mean,data_ROI,numberROIs,image_for_size)
%[D,Reg_all] = mean_time_course(data_for_mean, mask_ROIs,numberROIs)

size_data = size(data_mean);

[ROI_matrix] = resize_img_useTemp_imcalc(data_ROI,image_for_size);

vectorROI = reshape(ROI_matrix,size(ROI_matrix,1)*size(ROI_matrix,2)*size(ROI_matrix,3),1);

D = zeros(size_data(1,1),length(numberROIs));
Reg_all = zeros(length(numberROIs),1);
% mean time course extraction

for i = 1:size_data(1,1)
    try
        
    disp(['Extracting data for ' data_mean{i}]);
  
        X = spm_vol(data_mean{i});
        Y = spm_read_vols(X);
        
        if size(ROI_matrix)~=size(Y)
            [vol_matrix] = resize_img_useTemp_imcalc(data_mean{i},image_for_size);
            vector_vol = reshape(vol_matrix,size(vol_matrix,1)*size(vol_matrix,2)*size(vol_matrix,3),1);
        else
            vector_vol = reshape(Y,size(Y,1)*size(Y,2)*size(Y,3),1);
        end
        
        if numel(vectorROI)==numel(vector_vol)
            for j = 1:length(numberROIs)
                Reg_j = vector_vol(round(vectorROI) == numberROIs(j));
                Reg_all(j) = length(Reg_j);
                D(i,j)= mean(removenan_my(Reg_j,':'));
                %disp([num2str(length(Reg_j)) ' voxels found for file ' num2str(i)  ' region ID ' num2str(numberROIs(j)) ': ' num2str(D(i,j))]);
            end
        else
            disp(['Wrong dimensions: ' data_mean{i}]);
        end
    catch
        disp(['Failed: ' data_mean{i}]);
    end
end



