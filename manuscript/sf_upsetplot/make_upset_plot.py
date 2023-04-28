from upsetplot import from_memberships, plot
import matplotlib.pyplot as plt
import joblib

data = joblib.load("memberships.jl")
memb = from_memberships(data.values())
plot(memb, show_percentages=True, subset_size='sum', 
     min_subset_size=int(len(data.keys()) * 0.005),
     sort_by='cardinality')
plt.savefig("upsetplot.pdf")
