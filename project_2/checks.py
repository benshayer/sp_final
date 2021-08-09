import numpy.random as npr
import pandas as pd
import string
import random
import numpy as np

# nd1 = np.array([[2, 0, 1], [1, 1, 0], [7, 1, 1], [4, 0, 0]])
# nd2 = np.array([1, 2, 3, 4])

# ser1 = pd.Series(nd2, index=list('abcd'))
# ser2 = pd.Series(nd2, index=list('AbcD'))

# nd3=nd1.reshape(3,-4)
# print(nd3)

data = {'Country': ['Belgium', 'India', 'Brazil'], 'Capital': ['Brussels', 'New Delhi', 'Brasília'],
        'Population': [11190846, 1303171035, 207847528]}

df = pd.DataFrame(data,columns=['Country','Capital','Population'])
df=df.drop(1,axis=0)
print(df)
