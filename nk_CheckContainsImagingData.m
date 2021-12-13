function [contains_img, where_img] = nk_CheckContainsImagingData (NM, analind)

F = NM.analysis{analind}.params.TrainParam.FUSION.M;
nF = numel(F); where_img = false(1,nF);
for i=1:nF
    if strcmp(NM.datadescriptor{F(i)}.source, 'image')
        where_img(i) = true;
    end
end
contains_img=any(where_img);