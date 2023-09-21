# compute reference networks with mutual information as the similarity measure

import numpy as np
from sklearn.metrics import mutual_info_score
from pyspi.calculator import Calculator
from pyspi.data import Data

ref = np.asarray(ref)

refDat = Data(ref, dim_order = 'sp', normalise = False)

calc = Calculator(dataset=refDat, configfile = 'pyspiconfig.yaml')
calc.compute()
ref_conmat = calc.table['mi_gaussian'].to_numpy()

X = np.asarray(X)
n_rois = X.shape[1]
netwX = np.zeros((X.shape[0], int((n_rois**2-n_rois)/2)))

for i in range(X.shape[0]):

    iref = np.vstack([ref,X[i,:]])
    irefDat = Data(iref, dim_order = 'sp', normalise = True)

    calc = Calculator(dataset = irefDat, configfile = 'pyspiconfig.yaml')
    calc.compute()
    iref_conmat = calc.table['mi_gaussian'].to_numpy()

    dif_conmat = ref_conmat - iref_conmat

    netwX[i,:] = dif_conmat[np.triu_indices(n_rois, 1)]