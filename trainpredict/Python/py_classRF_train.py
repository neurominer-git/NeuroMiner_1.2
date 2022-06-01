# trainRFClassifier
# Python script to train a random forest classifier from NeuroMiner

from sklearn.ensemble import RandomForestClassifier
import pickle

rf = RandomForestClassifier(n_estimators = n_est, max_features = n_maxfeat, class_weight = clwght)
rf.fit(feat, lab)
model_file = f'{rootdir}/RFC_model.sav'
pickle.dump(rf, open(model_file, 'wb'))
