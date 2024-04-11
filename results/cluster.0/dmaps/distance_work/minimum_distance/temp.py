import pickle

with open("min_distances.pkl", "rb") as pklf:
    alist = pickle.load(pklf)
print(alist[0])
