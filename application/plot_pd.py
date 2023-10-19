import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

df = pd.read_csv("pd_manha.csv")
df.dtypes
df.loc[df["seqnames"]=="22","seqnames"] = 22
df.loc[df["seqnames"]=="X","seqnames"] = 23
df["seqnames"].astype("int")
print(df["seqnames"].unique())
df["lgABF"] = np.log10(df["ABF"])
df = df.sort_values(['seqnames','pos'])
df.reset_index(inplace=True, drop=True)
df["i"] = df.index
df.dtypes
plot = sns.relplot(data=df, x='i', y='lgABF', height = 6, aspect = 3.4,
                   hue='seqnames', palette = 'dark', legend=None,sizes=1) 
chrom_df=df.groupby('seqnames')['i'].median()
plot.ax.set_xticks(chrom_df)
plot.ax.set_xticklabels(chrom_df.index)
plot.ax.set_xlabel('CHR',fontsize="xx-large")
plot.ax.set_ylabel('lg ABF',fontsize="xx-large")
#plot.ax.set_ylim(0,20)
#plot.ax.axhline(y=20, linewidth = 2,linestyle="--",color="grey")
plot.fig.suptitle('Manhattan plot (pure)',fontsize="xx-large",weight='bold')
plt.savefig('pd.tiff')
plt.close("all")

df = pd.read_csv("comb_manha.csv")
df.dtypes
df.loc[df["seqnames"]=="22","seqnames"] = 22
df.loc[df["seqnames"]=="X","seqnames"] = 23
df["seqnames"].astype("int")
print(df["seqnames"].unique())
df["lgABF"] = np.log10(df["ABF"])
df = df.sort_values(['seqnames','pos'])
df.reset_index(inplace=True, drop=True)
df["i"] = df.index
df.dtypes
plot = sns.relplot(data=df, x='i', y='lgABF', height = 6, aspect = 3.4,
                   hue='seqnames', palette = 'dark', legend=None,sizes=1) 
chrom_df=df.groupby('seqnames')['i'].median()
plot.ax.set_xticks(chrom_df)
plot.ax.set_xticklabels(chrom_df.index)
plot.ax.set_xlabel('CHR',fontsize="xx-large")
plot.ax.set_ylabel('lg ABF',fontsize="xx-large")
#plot.ax.set_ylim(0,20)
#plot.ax.axhline(y=20, linewidth = 2,linestyle="--",color="grey")
plot.fig.suptitle('Manhattan plot (mixed)',fontsize="xx-large",weight='bold')
plt.savefig('comb.tiff')