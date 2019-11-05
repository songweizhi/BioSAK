import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
data = [1.5]*7 + [2.5]*2 + [3.5]*8 + [4.5]*3 + [5.5]*1 + [6.5]*8
print(data)
density = gaussian_kde(data)
xs = np.linspace(-2,10,200)
density.covariance_factor = lambda : .3
density._compute_covariance()
plt.plot(xs,density(xs))
plt.show()