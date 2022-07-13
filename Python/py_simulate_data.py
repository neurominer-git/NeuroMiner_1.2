# function that simulates data with the same structure as input data
import pandas as pd
import numpy as np
from sdv.tabular import GaussianCopula

#print(shape)
data = pd.read_csv(data_file)
#print(shape(data))
#labels = np.ndarray(labels)
data['label'] = labels
#model = GaussianCopula(primary_key = ID)
model = GaussianCopula()
#model = CopulaGAN()
model.fit(data)

sample = model.sample(n_obs)
#sample.head()
out_path = f'{rootdir}/simData.csv'
sample.to_csv(out_path,index=False)
