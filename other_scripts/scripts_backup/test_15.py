
import csv
import numpy

filename = '/Users/songweizhi/Desktop/test.csv'

# raw_data = open(filename)
#
#
# reader = csv.reader(raw_data, delimiter=',', quoting=csv.QUOTE_NONE)
# print(reader)
# x = list(reader)
# print(x)
# data = numpy.array(x).astype('float')
# print(data)
# print(data.shape)



# from numpy import loadtxt
# raw_data = open(filename, 'rt')
# data = loadtxt(raw_data, delimiter=",")
# print(data)
# print(data.shape)


from pandas import read_csv
from numpy import set_printoptions
from sklearn.preprocessing import MinMaxScaler

names = ['preg', 'plas', 'pres', 'skin', 'test', 'mass', 'pedi', 'age', 'class']
dataframe = read_csv(filename, names=names)
array = dataframe.values
X = array[:,:8]
Y = array[:,8]

scaler = MinMaxScaler(feature_range=(0, 1))
rescaledX = scaler.fit_transform(X)
set_printoptions(precision=3)

print(X[0:5,:])

print(rescaledX[0:5,:])




