import json
import seaborn as sb
import matplotlib.pyplot as plt

def plotter(fn, y=0, ax=None):
    full = json.load(open(fn))
    data = full['test']
    p = sb.histplot(ax=ax, data=data, x="perms", color='gray', edgecolor='gray', kde=False, stat='density')
    p = sb.kdeplot(ax=ax, data=data, x="perms", color='black')

    x = data['observed'] 
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
    if ax is None:
        ax = p
    ax.axvline(x, color='blue')
    ax.text(x, y, 'observed intersections', rotation=90, bbox=props, ma='center')
    p.set(xlabel="Intersection Count", ylabel="Permutation Density")
    print(data['alt'], round(data['p_val'], 3))
    return p

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))

random = 'novl'
p = plotter(f"results/promoters_perchrom{random}1000.json", 0.0011, ax[0, 0])
p.set(title="TR catalog intersection with promoters")

p = plotter(f"results/transcript_perchrom{random}1000.json", 1.25e-5, ax[0, 1])
p.set(title="TR catalog intersection with transcripts")

p = plotter(f"results/nonintron_perchrom{random}1000.json", 0.000075, ax[1, 0])
p.set(title="TR catalog intersection with non-introns")

p = plotter(f"results/intron_perchrom{random}1000.json", 1.5e-5, ax[1, 1])
p.set(title="TR catalog intersection with introns")

plt.tight_layout()
plt.savefig(f"TRGeneProm_{random}.pdf")
