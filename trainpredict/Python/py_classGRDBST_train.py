# train GRDBST Classifier
# Python script to train a gradient boosting classifier from NeuroMiner

from sklearn.ensemble import GradientBoostingClassifier
import pickle

gb = GradientBoostingClassifier(n_estimators = n_est,
                                loss = l,
                                learning_rate = lr,
                                subsample = subsamp,
                                max_depth = n_maxdepth)
gb.fit(feat, lab)
model_file = f'{rootdir}/GBC_model.sav'
pickle.dump(gb, open(model_file, 'wb'))
