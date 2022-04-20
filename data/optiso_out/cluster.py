import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt

# Data set
df = pd.read_csv('apc.tsv', sep='\t')
df = df.set_index('seq')
#df.drop(['fit'], axis=1, inplace=True)


# Default plot
sns.clustermap(df, metric="euclidean", standard_scale=1, method="ward")

# Show the graph
plt.show()
