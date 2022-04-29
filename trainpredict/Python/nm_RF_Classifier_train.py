# trainRFClassifier
# Python script to train a random forest classifier from NeuroMiner

from sklearn.ensemble import RandomForestClassifier
import pickle

def nm_RF_Classifier_train(train_features, train_labels, nest, nfeat):
    rf = RandomForestClassifier(n_estimators = nest, max_features = nfeat)
    rf.fit(train_features, train_labels)
    filename = 'RFC_model.sav'
    pickle.dump(rf, open(filename, 'wb'))
    return filename

model_file = nm_RF_Classifier_train(feat, lab, n_est, n_maxfeat)
