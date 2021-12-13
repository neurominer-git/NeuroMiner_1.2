function I = nk_RemapData(D, I, remap)

[~, m] = size(I.Y);
tY = nan(D.n_subjects_all, m);
tY(remap.vec2,:) = I.Y(remap.vec1,:);
t_ID = D.cases;
I.Y = tY;
I.ID = t_ID;
