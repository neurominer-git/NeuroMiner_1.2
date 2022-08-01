# trainRFRegressor
# Python script to train a random forest regressor from NeuroMiner

from sklearn.ensemble import RandomForestRegressor
import pickle

# criterion
if crit == 1:
    crit = 'squared_error' # if scikit-learn version <v1.0, then crit = 'mse'
elif crit == 2:
    crit = 'absolute_error'
elif crit == 3:
    crit = 'poisson'

# max_depth
if maxd == 0:
    maxd = None

# max_features
if n_maxfeat == -1:
    n_maxfeat = 'sqrt'
elif n_maxfeat == -2:
    n_maxfeat = 'log2'
elif n_maxfeat == 0:
    n_maxfeat = None

# min_samples_leaf
if minsl >= 1.0:
    minsl = int(minsl)

# min_samples_split
if minss > 1.0:
    minss = int(minss)

# max_leaf_nodes
if maxln == 0:
    maxln = None


# max_samples
if maxs == 0:
    maxs = None

rf = RandomForestRegressor(n_estimators = n_est,
        max_features = n_maxfeat,
        criterion = crit,
        max_depth = maxd,
        min_samples_split = minss,
        min_samples_leaf = minsl,
        min_weight_fraction_leaf = minwfl,
        max_leaf_nodes = maxln,
        min_impurity_decrease = minid,
        bootstrap = boot,
        oob_score = oobs,
        ccp_alpha = ccpa,
        max_samples = maxs
        ) # others: verbose, random_state, warm_start, n_jobs

rf.fit(feat, lab)
model_file = f'{rootdir}/RFR_model.sav'
pickle.dump(rf, open(model_file, 'wb'))
