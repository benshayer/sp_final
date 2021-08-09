from sklearn.cluster import KMeans
from sklearn.datasets import load_iris
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')
data_set = load_iris(return_X_y=True, as_frame=True)[0]
inertia_values = np.zeros(10)
for k in range(1, 11):
    inertia_values[k - 1] = KMeans(n_clusters=k, random_state=0).fit(data_set).inertia_
ki = [i for i in range(1, 11)]
plt.plot(ki, inertia_values, marker='.')
plt.xlabel("K values")
plt.ylabel("Inertia")
plt.annotate('elbow', xy=(2, inertia_values[1]),  xycoords='data',
            xytext=(0.5, 0.5), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),
            horizontalalignment='right', verticalalignment='top'
            )
plt.plot(range(1, 11), inertia_values, markevery=[1], marker='$◌$', markersize=30, label="circle", color="blue")
#plt.show()
plt.savefig('elbow.png')
