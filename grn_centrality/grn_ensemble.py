#Download dependencies before proceeding

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from arboreto.algo import grnboost2
import os

sns.set(rc={'figure.figsize':(12,8)})
sns.set(style='whitegrid')

def run_tf_search():
    wd = '/home/users/jvs15/MachineLearningProstate'
    #wd = '/Users/Jake/MachineLearningProstate'

    #Read in time-series RNA-seq data
    data_df = []
    names = []
    for name in os.listdir(f'{wd}/data/timeseries/'):
        data_df.append(pd.read_csv(f'{wd}/data/timeseries/{name}'))
        names.append(name)
    tf_df = pd.read_csv(f'{wd}/STAMPScreen/grn_centrality/Homo_sapiens_TF.csv')
    #List of TFs
    tfs = set(tf_df['Symbol'].tolist())
    #Create transcription factor DF 
    #Training can occur on just TFs or with all genes, depending on computational resources.
    onlytfs_df = []
    for i in data_df:
        x = i[i['Gene'].isin(tfs)]
        x = x.set_index('Gene')
        onlytfs_df.append(x)

    os.chdir(f'{wd}/STAMPScreen/grn_centrality')
    prdf_tot = []

    # Average out most important transcription factors

    for j in range(100):
        #Train network using GRNBoost2
        network = []
        ofile = []
        for i in range(len(onlytfs_df)):
            network.append(grnboost2(expression_data=onlytfs_df[i].T))
            print(f"Finished training {names[i]}")
            #Network file   
            ofile.append(f'trained_network_{names[i]}')
            network[i].to_csv(ofile[i])
            print("Saved network file!")

        #View Network
        network.clear()
        for i in range(len(ofile)):
            network.append(pd.read_csv(f'{wd}/STAMPScreen/grn_centrality/{ofile[i]}'))

        '''
        #Visualize edge weight distribution
        for i in range(len(network)):
            f, axes = plt.subplots(2, 1)

            sns.distplot(a=network[i]['importance'], norm_hist=True, ax=axes[0]).set(title= f'{names[i][0:len(names[i])-4]}')
            sns.distplot(a=network[i]['importance'].apply(np.log), norm_hist=True, ax=axes[1]).set(title= f'{names[i][0:len(names[i])-4]}')

            axes[0].set_xlabel('Connection weight')
            axes[1].set_xlabel('Connection weight (log)')
        '''

        #Create graph
        G = [] 
        for i in range(len(network)):
            G.append(nx.from_pandas_edgelist(df=network[i], source='TF', target='target', edge_attr='importance'))
            print('Loaded {:,} genes with {:,} edges.'.format(len(G[i].nodes), len(G[i].edges)))

        #Prune graph

        cutoff = 1

        print('Removing all edges with weight < {}...\n'.format(cutoff))

        for i in range(len(G)):
            bad_edges = [(s,t,w) for (s,t,w) in G[i].edges.data('importance') if w < cutoff]
            G[i].remove_edges_from(bad_edges)
            print('Graph now has {:,} genes and {:,} edges.'.format(len(G[i].nodes), len(G[i].edges)))

        # Run PageRank
        prdf = []
        for i in range(len(G)):
            pr = nx.pagerank(G[i], alpha=0.85, max_iter=50, weight='importance')
            #Create dataframe for PageRank values
            prdf.append(pd.DataFrame(pd.Series(pr)).reset_index())
            prdf[i].columns = ['Gene', 'PageRank']
            if len(prdf_tot) < len(names):
                prdf_tot.append(prdf[i])
            else:
                prdf_tot[i] = pd.concat([prdf_tot[i], prdf[i]])\
                    .groupby('Gene')['PageRank'].sum().reset_index()

        
        for i in range(len(prdf_tot)):
            prdf_tot[i].to_csv(f'aggregate_tfs_{names[i]}')
            print(prdf_tot[i].head())
            print(f'aggregate {names[i]} added')

        print(f"Finished round {j}")
    
    #View top TFs
    for i in range(len(prdf_tot)):
        sns.set(rc={'figure.figsize':(24,16)})
        sns.barplot(
            data=prdf_tot[i].sort_values('PageRank', ascending=False).head(45), 
            x='PageRank', 
            y='Gene'
        ).set(title=f'{names[i][0:len(names[i])-4]} 100 iterations')
        data=prdf_tot[i].sort_values('PageRank', ascending=False)
        #Save ranked TFs
        data.to_csv(f'aggregate_{names[i][0:len(names[i])-4]}_pageranked_tfs.csv')
        plt.show()




if __name__ == "__main__":
    run_tf_search()

