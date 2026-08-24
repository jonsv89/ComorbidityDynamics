import networkx as nx
import distanceclosure as dc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns


def create_weighted_network(fname, weight='weight', source='source', target='target'):
    '''
    fname (str): Name of the file containing the network as an edgelist table.
    weight (str): Column name of the attribute containing the relative risk edge weights
    source/target (str): Columns name containing the source/target disease code   
    '''

    gdf = pd.read_csv(fname, sep='\t')
    gdf = gdf[gdf[weight]>=1.01]
    gdf['distance'] = 1./gdf[weight]

    G = nx.from_pandas_edgelist(gdf, source=source, target=target, edge_attr=True, create_using=nx.DiGraph)        

    # Add color/category attribute
    cdf = pd.read_csv(f'disease_colors.txt', sep='\t')
    node_attr = dict(zip(cdf['disease'], cdf['category']))
    nx.set_node_attributes(G, values=node_attr, name="group")

    return G


def get_closure_mixing(G, node_attr="group"):
    '''
    G (nx.DiGraph): Networkx graph with distance edge weights.
    node_attr (str): Name of the node attribute to group.
    '''

    gc = dc.distance_closure(G, kind='ultrametric', weight='distance')
    mix_dict = nx.attribute_mixing_dict(gc, "group", normalized=False)

    G = nx.Graph()
    for u, edg in mix_dict.items():
        ku = np.sum(list(edg.values()))
        for v, count in edg.items():
            G.add_edge(u, v, count=count, no_count=ku-count)
    gdf = nx.to_pandas_edgelist(G)

    return gdf


def heatmap(exp_df, ref_df):

    or_df = pd.merge(exp_df, ref_df, on=['source', 'target'], how='inner', suffixes=('_exp', '_ref'))
    or_df['OR'] = (or_df['count_exp']*or_df['no_count_ref'])/(or_df['count_ref']*or_df['no_count_exp'])
    or_df['LogOR'] = or_df['OR'].apply(np.log10)

    or_df['SE'] = np.sqrt(or_df[['count_exp', 'count_ref', 'no_count_exp', 'no_count_ref']].rdiv(1).sum(axis=1))
    or_df['upper'] = or_df['LogOR']+1.96*or_df['SE']
    or_df['lower'] = or_df['LogOR']-1.96*or_df['SE']
    or_df['sign'] = np.sign(or_df['upper']*or_df['lower'])

    or_df = or_df[or_df['sign']==1.0]
    pivot_df = or_df.pivot(index='source', columns='target', values='LogOR')
    all_entities = sorted(list(set(pivot_df.index).union(set(pivot_df.columns))))

    return pivot_df.reindex(index=all_entities, columns=all_entities, fill_value=np.nan)


def plot_heatmap(heatmap_data, exp_name='Expected', ref_name='Reference', fig_title='Title'):

    fig, ax = plt.subplots(figsize=(16, 14))
        
    ax = sns.heatmap(heatmap_data, annot=True, cmap='PiYG', fmt='.2f', center=0, linewidths=0.8, linecolor='black')

    ax.set_xlabel('Prospective disease', fontsize=18)
    ax.set_ylabel('Index disease', fontsize=18)

    ax.tick_params(axis='y', labelrotation=0, labelsize=14)
    ax.tick_params(axis='x', labelrotation=270, labelsize=14)

    cbar = ax.collections[0].colorbar
    cbar.set_label('LogOdds Prospective Comorbidity', fontsize=16, labelpad=15)

    cbar.ax.set_title(exp_name, fontsize=16, pad=10)
    cbar.ax.set_xlabel(ref_name, fontsize=16, labelpad=10)

    plt.title(fig_title, fontsize=20)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    ## Time Window Comparison

    G1 = create_weighted_network('data/RR_net_both_0_1.txt', weight='RR_shrunk', source='disease_a', target='disease_b')
    G1c_df = get_closure_mixing(G1)

    G5 = create_weighted_network('data/RR_net_both_0_5.txt', weight='RR_shrunk', source='disease_a', target='disease_b')
    G5c_df = get_closure_mixing(G5)

    heatmap_data = heatmap(G1c_df, G5c_df)
    plot_heatmap(heatmap_data, '0-1 Window', '0-5 Window', 'Time Window Comparison')

    ## Catalunya Denmark Comparison

    D = create_weighted_network('data/denmark_raw.txt', weight='adjustedRR2.5', source='A', target='B')
    Dc_df = get_closure_mixing(D)

    heatmap_data = heatmap(G5c_df, Dc_df)
    plot_heatmap(heatmap_data, 'Catalunya', 'Denmark', 'Population Comparison')