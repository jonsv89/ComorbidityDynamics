"""
Propective Comorbidities Computation
"""


import networkx as nx
import distanceclosure as dc
import pandas as pd
import numpy as np

import argparse

import matplotlib.pyplot as plt
import seaborn as sns


def compute_closure(fname: str, attr_file: str,
                    weight: str='weight', source: str='source', target: str='target',
                    node: str='node', category: str='category',
                    save: bool=False) -> nx.DiGraph:
    '''
    Input:
        fname (str): Path to network saved as an edgelist.
        attr_file (str): Path to table of node category to be added to the network.
        
        weight (str, optional): Edgelist column label of `weight` attribute in `fname`. Default is 'weight'.
        source (str, optional): Edgelist column label of `source` attribute in `fname`. Default is 'source'.
        target (str, optional): Edgelist column label of `target` attribute in  `fname`. Default is 'target'.
        node (str, optional): Table column label of `node` attribute in `attr_file`. Default is 'node'.
        category (str, optional): Table column label of `category` attribute in `attr_file`. Default is 'category'.
        save (bool, optional): Whether to save intermediary files. Default is False.
        
    Output:
        G (nx.DiGraph): Networkx graph with distance edge weights.
    '''
       
    gdf: pd.DataFrame = pd.read_csv(filepath_or_buffer=fname, sep='\t')
    gdf = gdf[gdf[weight]>=1.01]
    gdf['distance'] = 1./gdf[weight]
    
    G: nx.DiGraph = nx.from_pandas_edgelist(df=gdf, source=source, target=target, edge_attr=True, create_using=nx.DiGraph)        

    # Add category attribute
    cdf: pd.DataFrame = pd.read_csv(filepath_or_buffer=attr_file, sep='\t')
    node_attr: dict = dict(zip(cdf[node], cdf[category]))
    nx.set_node_attributes(G, values=node_attr, name="group")
    
    G = dc.distance_closure(D=G, kind='ultrametric', weight='distance')
    
    if save:
        nx.write_edgelist(G, path=fname.replace('.txt', '_closure.txt'), data=True, delimiter='\t')   

    return G


def get_closure_mixing(G: nx.DiGraph, category: str='category') -> pd.DataFrame:
    
    '''
    Input:
        G (nx.DiGraph): Networkx graph with distance edge weights.
        category (str, optional): Table column label of `category` attribute in `attr_file`. Default is 'category'.

    Output:
        gdf (pd.DataFrame): Pandas pd.DataFrame containing the mixing matrix.
    '''

    mix_dict: dict = nx.attribute_mixing_dict(G, category, normalized=False)

    g: nx.Graph = nx.Graph()
    for u, edg in mix_dict.items():
        ku: float = np.sum(a=list(edg.values()))
        for v, count in edg.items():
            g.add_edge(u_of_edge=u, v_of_edge=v, count=count, no_count=ku-count)
    gdf: pd.DataFrame = nx.to_pandas_edgelist(g)
    
    return gdf


def pivot_table(pos_df: pd.DataFrame, neg_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Input:
        pos_df (pd.DataFrame): Pandas pd.DataFrame containing the mixing matrix of the positive network
        neg_df (pd.DataFrame): Pandas pd.DataFrame containing the mixing matrix of the negative network
    
    Output:
        pivot_df (pd.DataFrame): Pandas pd.DataFrame containing the pivoted heatmap data.
    '''

    or_df: pd.DataFrame = pd.merge(left=pos_df, right=neg_df, on=['source', 'target'], how='inner', suffixes=('_pos', '_neg'))
    or_df['OR'] = (or_df['count_pos']*or_df['no_count_neg'])/(or_df['count_neg']*or_df['no_count_pos'])
    or_df['LogOR'] = or_df['OR'].apply(func=np.log10)

    or_df['SE'] = np.sqrt(or_df[['count_pos', 'count_neg', 'no_count_pos', 'no_count_neg']].rdiv(1).sum(axis=1))
    or_df['upper'] = or_df['LogOR']+1.96*or_df['SE']
    or_df['lower'] = or_df['LogOR']-1.96*or_df['SE']
    or_df['sign'] = np.sign(or_df['upper']*or_df['lower'])

    or_df = or_df[or_df['sign']==1.0]
    pivot_df: pd.DataFrame = or_df.pivot(index='source', columns='target', values='LogOR')
    all_entities: list[str] = sorted(list(set(pivot_df.index).union(set(pivot_df.columns))))
    pivot_df = pivot_df.reindex(index=all_entities, columns=all_entities, fill_value=np.nan)

    return pivot_df


def plot_heatmap(heatmap_data: pd.DataFrame, pos_name: str='Positive', neg_name: str='Negative', fig_title: str='Title') -> None:
    '''
    Input:
        heatmap_data (pd.DataFrame): Pandas pd.DataFrame containing the pivoted heatmap data.
        pos_name (str, optional): Label for the positive network in the plot.
        neg_name (str, optional): Label for the negative network in the plot.
        fig_title (str, optional): Title for the plot.
    '''

    fig, ax = plt.subplots(figsize=(16, 14))
        
    ax = sns.heatmap(data=heatmap_data, annot=True, cmap='PiYG', fmt='.2f', center=0, linewidths=0.8, linecolor='black')

    ax.set_xlabel(xlabel='Prospective disease', fontsize=18)
    ax.set_ylabel(ylabel='Index disease', fontsize=18)

    ax.tick_params(axis='y', labelrotation=0, labelsize=14)
    ax.tick_params(axis='x', labelrotation=270, labelsize=14)

    cbar = ax.collections[0].colorbar
    cbar.set_label('LogOdds Prospective Comorbidity', fontsize=16, labelpad=15)

    cbar.ax.set_title(pos_name, fontsize=16, pad=10)
    cbar.ax.set_xlabel(neg_name, fontsize=16, labelpad=10)

    plt.title(label=fig_title, fontsize=20)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Process script inputs.")
    
    parser.add_argument("posfile", type=str, help="Path to network saved as an edgelist, whose edges will indicate positive logit.")
    parser.add_argument("negfile", type=str, help="Path to network saved as an edgelist, whose edges will indicate negative logit.")
    parser.add_argument("attrfile", type=str, help="Path to table of node category to be added to the network.")
    
    parser.add_argument("--posweight", type=str, default='weight', help="Edgelist column label of `weight` attribute in `posfile`.")
    parser.add_argument("--possource", type=str, default='source', help="Edgelist column label of `source` attribute in `posfile`.")
    parser.add_argument("--postarget", type=str, default='target', help="Edgelist column label of `target` attribute in `posfile`.")
    
    parser.add_argument("--negweight", type=str, default='weight', help="Edgelist column label of `weight` attribute in `negfile`.")
    parser.add_argument("--negsource", type=str, default='source', help="Edgelist column label of `source` attribute in `negfile`.")
    parser.add_argument("--negtarget", type=str, default='target', help="Edgelist column label of `target` attribute in `negfile`.")
    
    parser.add_argument("--attrnode", type=str, default='node', help="Table column label of `node` attribute in `attrfile`.")
    parser.add_argument("--catnode", type=str, default='category', help="Table column label of `category` attribute in `attrfile`.")
    
    parser.add_argument("--save", type=bool, default=False, help="Whether to save intermediary files.")
    parser.add_argument("--plot", type=bool, default=True, help="Whether to plot the results.")
    
    parser.add_argument("--pos_label", type=str, default='Positive Logit', help="Label for the positive network in the plot.")
    parser.add_argument("--neg_label", type=str, default='Negative Logit', help="Label for the negative network in the plot.")
    parser.add_argument("--fig_title", type=str, default='Prospective Comorbidities Comparison', help="Title for the plot.")

    args: argparse.Namespace = parser.parse_args()
    
    G1: nx.DiGraph = compute_closure(fname=args.posfile, attr_file=args.attrfile, 
                    weight=args.posweight, source=args.possource, target=args.postarget,
                    node=args.attrnode, category=args.catnode, save=args.save)
    
    G2: nx.DiGraph = compute_closure(fname=args.negfile, attr_file=args.attrfile, 
                        weight=args.negweight, source=args.negsource, target=args.negtarget,
                        node=args.attrnode, category=args.catnode, save=args.save)
    
    shared_nodes: set[str] = set(G1.nodes()).intersection(G2.nodes())
    G1 = G1.subgraph(nodes=shared_nodes).copy()
    G2 = G2.subgraph(nodes=shared_nodes).copy()
    
    comparison_data: pd.DataFrame = pivot_table(pos_df=get_closure_mixing(G=G1, category=args.catnode),
                                             neg_df=get_closure_mixing(G=G2, category=args.catnode))
    if args.save:
        comparison_data.to_csv(path_or_buf='comparison_pivot_table.txt', sep='\t', index=True) 
    
    if args.plot:
        plot_heatmap(heatmap_data=comparison_data, pos_name=args.pos_label, neg_name=args.neg_label, fig_title=args.fig_title)
    else:
        print(comparison_data)