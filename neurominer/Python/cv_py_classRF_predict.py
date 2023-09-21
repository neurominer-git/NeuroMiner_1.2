from sklearn.ensemble import RandomForestClassifier
import pickle
from scipy.io import savemat
import os

rf = pickle.load(open(model_name, 'rb'))
#os.remove(model_name)
predictions = rf.predict(test_feat)
probabilities = rf.predict_proba(test_feat)

res_dict = {'predictions': predictions, 'probabilities': probabilities}

py_modeldir = os.path.join(rootdir,"Py_modelfiles")
if not os.path.isdir(py_modeldir):
    os.makedirs(py_modeldir)
    print(f"created folder: {py_modeldir}") 

results_file = os.path.join(py_modeldir,"RFpredict_output.mat")

savemat(results_file, res_dict)
