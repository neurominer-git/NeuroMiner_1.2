function [D,Reg_all] = mean_time_course_NM(data_mean,data_ROI,numberROIs,image_for_size,dir_save)
%[D,Reg_all] = mean_time_course(data_for_mean, mask_ROIs,numberROIs)
global JSMEM
size_data = size(data_mean);


% resize atlas to image_for_size (NM.brainmask)
if ~isfield(JSMEM, 'ROI_matrix') || isempty(JSMEM.ROI_matrix)
    [ROI_matrix] = resize_img_useTemp_imcalc_NM(data_ROI,image_for_size,dir_save);
    JSMEM.ROI_matrix = ROI_matrix;
else
    [ROI_matrix] = JSMEM.ROI_matrix;
end
vectorROI = reshape(ROI_matrix,size(ROI_matrix,1)*size(ROI_matrix,2)*size(ROI_matrix,3),1);

D = zeros(size_data(1,1),length(numberROIs));
Reg_all = zeros(length(numberROIs),1);
% mean time course extraction

for i =1:size_data(1,1)
    try
        
    disp(['Extracting data for ' data_mean{i}]);
  
        X = spm_vol(data_mean{i});
        Y = spm_read_vols(X);
        
        if size(ROI_matrix)~=size(Y)
            [ROI_matrix] = resize_img_useTemp_imcalc_NM(data_ROI,data_mean{i},dir_save);
            vectorROI = reshape(ROI_matrix,size(ROI_matrix,1)*size(ROI_matrix,2)*size(ROI_matrix,3),1);
        end
        vector_vol = reshape(Y,size(Y,1)*size(Y,2)*size(Y,3),1);
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



%     for i = 1:size_data(1,1)
%         try
%             % transpose data_mean row to vector column
%             % Y = data_mean(i,:)';
%             X = spm_vol(data_mean{i});
%             Y = spm_read_vols(X);
%             
%             %         if size(vectorROI)~=size(Y)
%             %             disp("Atlas and image don't have the same voxel size")
%             %         end
%             %         vector_vol = Y;
%             
%             if size(ROI_matrix)~=size(Y)
%                 [ROI_matrix] = resize_img_useTemp_imcalc(data_ROI,data_mean{i});
%                 vectorROI = reshape(ROI_matrix,size(ROI_matrix,1)*size(ROI_matrix,2)*size(ROI_matrix,3),1);
%             end
%             vector_vol = reshape(Y,size(Y,1)*size(Y,2)*size(Y,3),1);
%             
%             
%             if numel(vectorROI)==numel(vector_vol)
%                 for j = 1:length(numberROIs)
%                     Reg_j = vector_vol(round(vectorROI) == numberROIs(j));
%                     Reg_all(j) = length(Reg_j);
%                     D(i,j)= mean(removenan_my(Reg_j,':'));
%                     %disp([num2str(length(Reg_j)) ' voxels found for file ' num2str(i)  ' region ID ' num2str(numberROIs(j)) ': ' num2str(D(i,j))]);
%                 end
%             else
%                 disp(['Wrong dimensions: ' data_mean(i)]);
%             end
%         catch
%             disp(['Failed: ' data_mean(i)]);
%         end
%     end
% else
%     D = JSMEM.data_PET;
end


