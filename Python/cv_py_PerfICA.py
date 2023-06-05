from sklearn.decomposition import FastICA
import numpy as np
import pickle
import uuid

# train 
if mode == "train":
    if num_ics == 0:
        ica = FastICA(whiten="arbitrary-variance", random_state = 1234)
    else: 
        ica = FastICA(n_components=num_ics, whiten="arbitrary-variance", random_state = 1234)
    ica_model = ica.fit(data)  # Reconstruct signals
    S = ica_model.transform(data)
    ICs = ica_model.components_
    random_name = uuid.uuid4().hex
    model_file = f"{rootdir}/ica_model_{random_name}.sav"
    pickle.dump(ica, open(model_file, "wb"))
# test
elif mode == "test": 
    ica = pickle.load(open(ica_model, "rb"))
    S = ica.transform(data)
elif mode == "inverse_transform":
    data = np.asarray(data)
    print(data.ndim)
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    ica = pickle.load(open(ica_model, "rb"))
    S = ica.inverse_transform(data)
                      





