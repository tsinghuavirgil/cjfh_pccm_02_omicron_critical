import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm, rcParams
##########
prm_mtx = pd.read_csv('pseudo_gene.csv',header = 0, index_col=0)
sns.set_style('white')
sns.lmplot(scatter=False,data=prm_mtx,x='pseudotime',y='expression',
           hue='genes',
           #palette = neu_subset_palette,
           order=2,        
          )
plt.savefig('prm_gene.pdf')