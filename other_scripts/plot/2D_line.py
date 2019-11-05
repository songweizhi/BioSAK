import matplotlib.pyplot as plt
year = [1950, 1970, 1990, 2010]
pop = [2.5,3.7, 5.2, 7.0]
pop_2 = [1.5,2.7, 3.2, 6.0]
plt.plot(year,pop)
plt.scatter(year,pop)
plt.plot(year,pop_2)
plt.scatter(year,pop_2)
plt.show()

