# trainRFRegressor
# Python script to train a random forest regressor from NeuroMiner

from sklearn.ensemble import RandomForestRegressor
import pickle

rf = RandomForestRegressor(n_estimators = n_est, max_features = n_maxfeat)
rf.fit(feat, lab)
model_file = f'{rootdir}/RFC_model.sav'
pickle.dump(rf, open(model_file, 'wb'))
