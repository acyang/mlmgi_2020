# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 04:07:57 2020

@author: 1203087
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Importing the dataset
dataset = pd.read_csv('data.csv')
X = dataset.iloc[:,1:2].values
y = dataset.iloc[:,2].values

# Splitting the dataset into the Training set and Test set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)

# Fitting the Regression Model to the dataset
from sklearn.ensemble import RandomForestRegressor
regressor = RandomForestRegressor(n_estimators=100, random_state=0, )
regressor.fit(X_train, y_train)

# Visualising the Regression results for higher and smoother curve
X_grid = np.arange(min(X_train), max(X_train), 0.01)
X_grid = X_grid.reshape(len(X_grid), 1)
plt.scatter(X_train, y_train, color = 'red')
plt.plot(X_grid, regressor.predict(X_grid), color='blue')
plt.title('Random Forest Regression')
plt.xlabel('Density [eV]')
plt.ylabel('Band Gap [eV]')

#density of MoS2 : 5.06 g/cm³
test = np.array([5.06]).reshape(len(test), 1)
print("The caculated bandgap is: 1.8eV", )
print("The predictive bandgap is: ",regressor.predict(test))
plt.scatter(test, regressor.predict(test), marker='^', color='green')
plt.show()
