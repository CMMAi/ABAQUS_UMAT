import pandas as pd
import matplotlib.pyplot as plt

### Read data ###
dataname = '001mm_bottom'
data = pd.read_csv(dataname+'.csv')

X = data['time'].values
colname = ['E11','E12','E22','S11','S12','S22']
for i in range(3):
    strain = colname[i]+'_b'
    stress = colname[i+3]+'_b'

    ### plot ###
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10,8))
    plt.plot(data[strain].values,data[stress].values)
    #plt.scatter(X, Y, c='b', marker='o', s=15)
    plt.title(colname[i]+' & '+colname[i+3])
    plt.xlabel('Strain', fontsize=20)
    plt.ylabel('Stress', fontsize=20)
    plt.savefig(colname[i]+' & '+colname[i+3])

    ### plot ###
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10,8))
    plt.plot(X,data[strain].values)
    plt.title(colname[i])
    plt.xlabel('time', fontsize=20)
    plt.ylabel('Strain', fontsize=20)
    plt.savefig(colname[i])

    ### plot ###
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10,8))
    plt.plot(X,data[stress].values)
    plt.title(colname[i+3])
    plt.xlabel('time', fontsize=20)
    plt.ylabel('Stress', fontsize=20)
    plt.savefig(colname[i+3])
