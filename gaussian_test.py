from imports import *

df = pd.read_csv("./.gaussian_spread_data/results.diff",delimiter=" ")

mags = df["Variable"].to_numpy()

mean, stdev = sp.stats.norm.fit(mags)

plt.figure(figsize=(8,6), dpi=200)
n = plt.hist(mags, bins=100)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 1000)
y = sp.stats.norm.pdf(x, mean, stdev)
plt.plot(x, y)
plt.xlabel("Apparent magnitude of TYC 4502-724-1")
plt.ylabel("Counts per bin")
plt.tight_layout()
plt.savefig("gaussian_long_axis.jpeg")
plt.show()

plt.figure(figsize=(8,6), dpi=200)
plt.xlim([xmin,(mean + (mean-xmin))])
plt.hist(mags, bins=100,color='b')
plt.plot(x, y,c='r',linestyle='dashed')
plt.xlabel("Apparent magnitude of TYC 4502-724-1")
plt.ylabel("Counts per bin")
plt.tight_layout()
plt.savefig("gaussian_short_axis.jpeg")
plt.show()

print("Mean Magnitude +/- Standard Deviation: ")