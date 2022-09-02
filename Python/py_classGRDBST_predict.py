'''
Author: Clara Vetter 
Last changed: 02.09.2022

This script is called by nk_GetTestPerf_GRDBST.m and computes the test 
performance of a gradient boosting classification model on the test sample. 

Input from MATLAB: 
    - model_name = path to file containing the model's information from 
    training
    - test_feat = the test set (test X)

Output: 
    - the predictions are stored in a MATLAB file GBpredict_output.mat in 
    the analysis' folder (which will be loaded in by Matlab to proceed) 
'''

from sklearn.ensemble import GradientBoostingClassifier
import pickle
from scipy.io import savemat
import uuid
#import os

gb = pickle.load(open(model_name, 'rb'))
#os.remove(model_name)
predictions = gb.predict(test_feat)
probabilities = gb.predict_proba(test_feat)

res_dict = {'predictions': predictions, 'probabilities': probabilities}

results_file = f'{rootdir}/GBpredict_output.mat'

savemat(results_file, res_dict)
