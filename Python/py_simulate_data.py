# function that simulates data with the same structure as input data
import pandas as pd
import numpy as np
from sdv.tabular import GaussianCopula

data = pd.read_csv(data_file)
data['label'] = labels

# if n_obs is an array with more than 1 element, each element specifies the amount of observations to simulate within each group (i.e. same label)
if type(n_obs) is not int:
    # check first whether len(n_obs) matches the number of unique labels
    if len(n_obs) != len(np.unique(labels)):
        raise ValueError('Number of observations to simulate does not match number of groups')

    # preallocate space for simulated dataset
    sim_sample = pd.DataFrame(columns = data.columns)
    for i in range(len(np.unique(labels))):
        # split dataset
        group_df = data[data['label'] == np.unique(labels)[i]]
        model = GaussianCopula()
        model.fit(group_df)
        sim_group = model.sample(n_obs[i])
        sim_sample = sim_sample.append(sim_group)
else:
    # model = GaussianCopula(primary_key = ID)
    model = GaussianCopula()
    # model = CopulaGAN()
    model.fit(data)
    sim_sample = model.sample(n_obs)

out_path = f'{rootdir}/simData.csv'
sim_sample.to_csv(out_path,index=False)
