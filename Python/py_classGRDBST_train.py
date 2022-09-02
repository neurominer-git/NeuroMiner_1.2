'''
Author: Clara Vetter 
Last changed: 02.09.2022

This script is called by nk_GetParam2_GRDBST.m and trains a gradient 
boosting classification model on the training data. 

Input from MATLAB: 
    - parameters: 
        n_est
        l
        lr
        subsamp
        n_maxdepth
    - feat = training data
    - lab = label 
    - rootdir = path to analysis directory

Output: 
    - the model is stored in the analysis directory in a .sav file that 
    will later be read in by the corresponding test script 
    (py_classGRDBST_predict.py). The file name is unique.
    
'''


from sklearn.ensemble import GradientBoostingClassifier
import pickle
import uuid

gb = GradientBoostingClassifier(n_estimators = n_est,
                                loss = l,
                                learning_rate = lr,
                                subsample = subsamp,
                                max_depth = n_maxdepth)
gb.fit(feat, lab)
#model_file = f'{rootdir}/GBC_model_{n_est}_{l}_{lr}_{subsamp}_{n_maxdepth}.sav'
random_name = uuid.uuid4().hex;
model_file = f'{rootdir}/GBC_model_{random_name}.sav';
pickle.dump(gb, open(model_file, 'wb'))
