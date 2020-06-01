# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/01/2020"
__description__ = "This script compute rank TOPSIS from CMIP5 models end OBS basedata"

import pandas as pd
import numpy as np

# Import data Prec and Temp(AMZ - NEB - MATOPIBA)
data = pd.read_csv('temp_amz_cmip5.csv')
data = data.values[:,1:]

w=np.array([0.25, 0.25, 0.25, 0.25])
i=np.array([1, 1, 1, 1])

# Step 1 -  Normalize data
# Taking axis=0 to sum values along column
normalizationFactor = np.sqrt(np.sum(data**2, axis=0, dtype=float), dtype=float)
#~ print(normalizationFactor.shape)

# Broadcasting operation to divide Xij with normalization factor
normalizedData = (data/normalizationFactor)

#Rounding normalized data values to 3 decimal places
normalizedData = np.round(normalizedData.astype(np.float64),decimals=2)
#~ print("Normalized Data:", normalizedData)

# Step  2 - Multiple each evaluation by the associated weigth:
wgtNormalizedData = normalizedData*w
#~ print(w.shape)
#~ print(wgtNormalizedData)

idealBest = []
idealWorst = []

# Step3 - Positive and negative idea solution
for x in range(data.shape[1]):
	if i[x]==1:
		idealBest.append(max(wgtNormalizedData[:,x]))
		idealWorst.append(min(wgtNormalizedData[:,x]))
	if i[x]==0:
		idealBest.append(min(wgtNormalizedData[:,x]))
		idealWorst.append(max(wgtNormalizedData[:,x]))

# Step 4 - Determine the distance to the negative and positive ideal solution 
distanceFromBest = np.sqrt(np.sum((wgtNormalizedData-idealBest)**2, axis=1, dtype=float), dtype=float)
distanceFromBest = distanceFromBest.reshape(distanceFromBest.shape[0], -1)
#~ print("DB",distanceFromBest)

distanceFromWorst = np.sqrt(np.sum((wgtNormalizedData-idealWorst)**2, axis=1, dtype=float), dtype=float)
distanceFromWorst = distanceFromWorst.reshape(distanceFromWorst.shape[0], -1)
#~ print("DW",distanceFromWorst)

# Step 5 - Calculate the relative closeness to the ideal solution
totalDistance = distanceFromBest+distanceFromWorst
#~ print("TD",totalDistance)

performance = distanceFromWorst/totalDistance
print((performance.tolist()).sort)

order = performance.argsort(axis=0)
print(order)

ranks = order.argsort(axis=0)
#~ print(ranks)

# Converting ranks to 1-d numpy array
ranks=ranks.reshape(ranks.shape[0],)

# Print rank table
print('Item', 'Rank', sep='\t')
for idx,x in enumerate(ranks):
	print(idx+1, ranks.shape[0]-(x), sep='\t', end='\n')
