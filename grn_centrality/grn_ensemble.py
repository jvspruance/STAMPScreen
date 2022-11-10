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
    tf_list = []
    for name in os.listdir(f'{wd}/data/timeseries/'):
        # Read in csv from expression data and add to list
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
        tf_list.append(list(x.index))

    os.chdir(f'{wd}/STAMPScreen/grn_centrality')
    prdf_tot = []

    ### Definitely needs edits--make dataframes all for one dose, not for individual
    newtfdf = []
    for i in range(len(data_df[0].columns)):
        newtfdf.append(pd.concat([x.iloc[:, i] for x in data_df], axis=1))
    doses = [x.split('_')[2] if (len(x.split('_')) > 2) else x.split('_')[1] for x in rows]




    # Average out most important transcription factors

    for j in range(1):
        #Train network using GRNBoost2
        network = []
        ofile = []
        for i in range(len(onlytfs_df)):
            network.append(grnboost2(expression_data=onlytfs_df[i].T))
            # x = data_df[i]
            # x = x.set_index('Gene')
            # network.append(grnboost2(expression_data=x.T, tf_names = tf_list[i]))
            print(f"Finished training {names[i]}")
            #Network file   
            ofile.append(f'trained_network_{names[i]}')
            network[i].to_csv(ofile[i])
            print("Saved network file!")

        #View Network
        network.clear()
        for i in range(len(ofile)):
            network.append(pd.read_csv(f'{wd}/STAMPScreen/grn_centrality/{ofile[i]}'))

        #Create graph
        G = make_graphs(network)
        # for i in range(len(network)):
        #     G.append(nx.from_pandas_edgelist(df=network[i], source='TF', target='target', edge_attr='importance'))
        #     print('Loaded {:,} genes with {:,} edges.'.format(len(G[i].nodes), len(G[i].edges)))

        #Prune graph

        cutoff = 1
        G = prune_graphs(G, cutoff)
        # print('Removing all edges with weight < {}...\n'.format(cutoff))

        # for i in range(len(G)):
        #     bad_edges = [(s,t,w) for (s,t,w) in G[i].edges.data('importance') if w < cutoff]
        #     G[i].remove_edges_from(bad_edges)
        #     print('Graph now has {:,} genes and {:,} edges.'.format(len(G[i].nodes), len(G[i].edges)))

        # Run PageRank -- keep track of aggregated data in prdf_tot
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

        # Store aggregated pagerank data 
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



def make_graphs(networks):
    """
    networks: list of networks
    returns a list of graphs
    """
    G = [] 
    for i in range(len(networks)):
        G.append(nx.from_pandas_edgelist(df=networks[i], source='TF', target='target', edge_attr='importance'))
        print('Loaded {:,} genes with {:,} edges.'.format(len(G[i].nodes), len(G[i].edges)))
    return G

def prune_graphs(graphs, cutoff):
    """
    graphs: list of graphs
    cutoff: int. Edge weights greater than cutoff are kept 
    returns pruned graphs
    """
    print('Removing all edges with weight < {}...\n'.format(cutoff))

    for i in range(len(graphs)):
        bad_edges = [(s,t,w) for (s,t,w) in graphs[i].edges.data('importance') if w < cutoff]
        graphs[i].remove_edges_from(bad_edges)
        print('Graph now has {:,} genes and {:,} edges.'.format(len(graphs[i].nodes), len(graphs[i].edges)))
    return graphs

def run_pagerank(G):
    """
    G: list of graphs (ideally pruned)
    returns list of pagerank dataframes
    """



if __name__ == "__main__":
    run_tf_search()

