import pandas as pd
import numpy as np

print("pd")
df = pd.read_csv("pd_manha.csv")
df.loc[df["seqnames"]=="22","seqnames"] = 22
df.loc[df["seqnames"]=="X","seqnames"] = 23
df["seqnames"].astype("int")
df = df.sort_values("ABF",ascending=False)
group = df.groupby("seqnames")
q = pd.DataFrame()
for i in range(1,24):
    print(i)
    p = group.get_group(i)
    q = q.append(p.head(20))
    print(q.head())
q.to_csv("sum_pd.csv",index=False)

print("comb")
df = pd.read_csv("comb_manha.csv")
df.loc[df["seqnames"]=="22","seqnames"] = 22
df.loc[df["seqnames"]=="X","seqnames"] = 23
df["seqnames"].astype("int")
df = df.sort_values("ABF",ascending=False)
group = df.groupby("seqnames")
q = pd.DataFrame()
for i in range(1,24):
    print(i)
    p = group.get_group(i)
    q = q.append(p.head(20))
    print(q.head())
q.to_csv("sum_comb.csv",index=False)