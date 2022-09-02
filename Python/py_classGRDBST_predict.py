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
