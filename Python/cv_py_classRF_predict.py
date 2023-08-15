from sklearn.ensemble import RandomForestClassifier
import pickle
from scipy.io import savemat
#import os

rf = pickle.load(open(model_name, 'rb'))
#os.remove(model_name)
predictions = rf.predict(test_feat)
probabilities = rf.predict_proba(test_feat)

res_dict = {'predictions': predictions, 'probabilities': probabilities}

if not os.path.isdir(f'{rootdir}/Py_modelfiles')
    os.makedirs(f'{rootdir}/Py_filesmodel')
    print(f'created folder: {rootdir}/Py_modelfiles') 

results_file = f'{rootdir}/Py_modelfiles/RFpredict_output.mat'

savemat(results_file, res_dict)
