from sklearn.decomposition import FastICA
import numpy as np
import pickle
import uuid

# train 
if mode == 'train':
    ica = FastICA(n_components=3, whiten="arbitrary-variance")
    S = ica.fit_transform(data)  # Reconstruct signals
 
    random_name = uuid.uuid4().hex;
    model_file = f'{rootdir}/ica_model_{random_name}.sav';
    pickle.dump(ica, open(model_file, 'wb'))
# test
else: 
    ica = pickle.load(open(ica_model, 'rb'))
    S = ica.inverse_transform(data)





