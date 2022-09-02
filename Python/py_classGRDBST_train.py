# train GRDBST Classifier
# Python script to train a gradient boosting classifier from NeuroMiner

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
