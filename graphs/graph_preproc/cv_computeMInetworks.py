# compute reference networks with mutual information as the similarity measure

import numpy as np
from sklearn.metrics import mutual_info_score
from pyspi.calculator import Calculator

calc = Calculator(dataset=ref, configfile = 'pyspiconfig.yaml')
calc.compute()
ref_conmat = calc.table['mi_gaussian'].to_numpy()

n_rois = X.shape[1]
netwX = np.zeros((X.shape[0], int((n_rois**2-n_rois)/2)))

for i in range(X.shape[0]):

    iref = pd.concat([ref,X.iloc[[i]]])

    calc = Calculator(dataset = iref, configfile = 'pyspiconfig.yaml')
    calc.compute()
    iref_conmat = calc.table['mi_gaussian'].to_numpy()

    dif_conmat = ref_conmat - iref_conmat

    netwX[i,:] = dif_conmat[np.triu_indices(n_rois, 1)]
