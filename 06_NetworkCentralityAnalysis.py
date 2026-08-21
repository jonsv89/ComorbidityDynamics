#!/usr/bin/env python3
import argparse

import networkx as nx
import pandas as pd
import numpy as np
import os


def compute_metrics(file, pathLoad = '', pathSave = '', diseasesOut = ['M54', 'J00', 'T14'], edgeWeight = ['RR_shrunk', 'cases_event', 'RR_shrunk x cases_event'], verbose = True):
    """
    Compute network metrics (weighted outdegree and pagerank) using different edge-weighting strategies.

    This function generates and evaluates disease networks from the input dataset using one or more edge-weight definitions. 
    Specified diseases can be excluded from the analysis prior to network construction.
    Results are saved to disk for downstream analyses.
    
    Parameters
    ----------
    file : str
        Input file containing the disease association network data.

    pathLoad : str, optional, default=''
        Directory containing the input file.
        
    pathSave : str, optional, default=''
        Directory where to save the results.

    diseasesOut : list of str, optional
        List of disease codes to exclude from the network.
        Default is ['M54', 'J00', 'T14'].

    edgeWeight : list of str, optional
        Edge weighting schemes to evaluate.
        Default is ['RR_shrunk', 'cases_event', 'RR_shrunk x cases_event'].

    verbose : bool, optional, default=False
        If True, display progress messages and additional information during processing.

    Save
    -------
    File:
        The computed network metrics are saved in the path of pathSave with the same file name as the input name.
    """

    #Checking filepahts are correct
    filepathLoad = os.path.join(pathLoad, file)
    if not os.path.exists(filepathLoad):
        raise FileNotFoundError(f"'{filepathLoad}' does not exist. Please provide a valid input file.")

    if not os.path.exists(pathSave):
        raise FileNotFoundError(f"'{pathSave}' does not exist. Please provide a valid path to save the results.")
        

    #Load the data 
    ext = file.split('.')
    if ext == ".csv": sep = ","
    elif ext in [".txt", ".tsv"]: sep = "\t"
    df = pd.read_csv(filepathLoad, sep=sep)
    
    #TakingOut the diseases
    if verbose: 
        diseasesA = df['disease_a'].unique()
        diseasesB = df['disease_b'].unique()
        diseases = np.unique(np.concatenate((diseasesA, diseasesB)))
        print('Diseases in A: {}, B: {}, both: {}'.format(len(diseasesA), len(diseasesB), len(diseases)))
    if len(diseasesOut):
        df = df[~df['disease_b'].isin(diseasesOut)]
    diseasesA = df['disease_a'].unique()
    diseasesB = df['disease_b'].unique()
    diseases = np.unique(np.concatenate((diseasesA, diseasesB)))
    if verbose: print('Diseases before removing {} in A: {}, B: {}, both: {}'.format(diseasesOut, len(diseasesA), len(diseasesB), len(diseases)))


    dfRes = pd.DataFrame(index= diseases)
    for weightFeature in edgeWeight:
        if verbose: print('Starting computing of {}'.format(weightFeature))
            
        if weightFeature == ' x '.join([edgeWeight[0], edgeWeight[1]]):
            df.insert(2, weightFeature, [item[weightFeature.split(' x ')[0]]*item[weightFeature.split(' x ')[1]] for ind, item in df.iterrows()])
            
        if weightFeature in df.columns: 
            #Creating the network
            G = nx.from_pandas_edgelist(df[['disease_a', 'disease_b', weightFeature]], 
                                        source='disease_a', target='disease_b', 
                                        edge_attr=weightFeature,
                                        create_using=nx.DiGraph())

            #Calculating stats
            outdegree_weight = dict(G.out_degree(weight=weightFeature)) #outdegree with the weight
            pr = nx.pagerank(G, weight=weightFeature) #pagerank

            #Saving the results
            dfToSaveODW = pd.DataFrame([outdegree_weight], index=[weightFeature+'_OutDegreeW']).transpose()
            dfRes = pd.merge(dfRes, dfToSaveODW, how='left', left_index=True, right_index=True)

            dfToSavePG = pd.DataFrame([pr], index=[weightFeature+'_PageRank']).transpose()
            dfRes = pd.merge(dfRes, dfToSavePG, how='left', left_index=True, right_index=True)
        else:
            print('{} not in the input file'.format(weightFeature))

    dfRes = dfRes.reset_index().drop_duplicates().set_index('index')
    if verbose: print('Computing Done')

    saveFileName = file.split('.')[0]
    dfRes.to_csv(pathSave+'{}_metrics.csv'.format(saveFileName))
    
    
def compute_rankings(file, pathLoad = '', pathSave = '', metrics = ['OutDegreeW', 'PageRank'], edgeWeight = ['RR_shrunk', 'cases_event', 'RR_shrunk x cases_event'], verbose = True):
    """
    Compute node rankings in a disease network using multiple network centrality metrics and edge-weighting schemes.
 
    Calculates disease rankings based on the selected network metrics. 
    Rankings are computed independently for each edge-weight definition, allowing comparison of how
    different weighting strategies influence node importance. 
    Results are saved to disk for downstream analyses.

    Parameters
    ----------
    file : str
        Input file containing the disease association network data.

    pathLoad : str, optional, default=''
        Directory containing the input file.
        
    pathSave : str, optional, default=''
        Directory where to save the results.

    metrics : list of str, optional
        Network metrics used to rank nodes.
        Default is ['OutDegreeW', 'PageRank'].

    edgeWeight : list of str, optional
        Edge weighting schemes to evaluate.
        Default is ['RR_shrunk', 'cases_event', 'RR_shrunk x cases_event'].

    verbose : bool, optional, default=False
        If True, display progress messages and additional information during processing.

    Save
    -------
    File:
        The computed network metrics are saved in the path of pathSave with the same file name as the input name.
    """
    
    #Checking filepahts are correct
    filepathLoad = os.path.join(pathLoad, file)
    if not os.path.exists(filepathLoad):
        raise FileNotFoundError(f"'{filepathLoad}' does not exist. Please provide a valid input file.")
    
    if not os.path.exists(pathSave):
        raise FileNotFoundError(f"'{pathSave}' does not exist. Please provide a valid path to save the results.")
        
    #Load the data 
    df = pd.read_csv(filepathLoad, index_col='index')
    
    dfRanking = pd.DataFrame(index=df.index)
    for i, weightFeature in enumerate(edgeWeight):

        for metric in metrics:
            col = weightFeature +'_'+ metric
            
            if col in df.columns:
                if verbose: print('Starting computing of {}'.format(col))
                dfTmp = df[col].sort_values(ascending=False).reset_index().reset_index().set_index('index')[['level_0']]
                dfTmp.columns = [weightFeature+'_'+metric]
                dfRanking = pd.merge(dfRanking, dfTmp, right_index=True, left_index=True)
            else:
                print('{} not in the input file'.format(col))
    
    if verbose: print('Computing Done')                
                     
    saveFileName = '_'.join(file.split('.')[0].split('_')[:-1])
    dfRanking.to_csv(pathSave+'{}_ranking.csv'.format(saveFileName))
    

def main():

    parser = argparse.ArgumentParser(
        description="Disease network analysis"
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    # -------------------------
    # compute_metrics
    # -------------------------

    parser_metrics = subparsers.add_parser(
        "metrics",
        help="Compute network metrics"
    )

    parser_metrics.add_argument("file")

    parser_metrics.add_argument(
        "--pathLoad",
        default=""
    )

    parser_metrics.add_argument(
        "--pathSave",
        default=""
    )

    parser_metrics.add_argument(
        "--diseasesOut",
        nargs="+",
        default=['M54', 'J00', 'T14']
    )

    parser_metrics.add_argument(
        "--edgeWeight",
        nargs="+",
        default=['RR_shrunk', 'cases_event', 'RR_shrunk x cases_event']
    )

    parser_metrics.add_argument(
        "--verbose",
        action="store_true"
    )

    # -------------------------
    # compute_rankings
    # -------------------------

    parser_rankings = subparsers.add_parser(
        "rankings",
        help="Compute node rankings"
    )

    parser_rankings.add_argument("file")

    parser_rankings.add_argument(
        "--pathLoad",
        default=""
    )

    parser_rankings.add_argument(
        "--pathSave",
        default=""
    )

    parser_rankings.add_argument(
        "--metrics",
        nargs="+",
        default=['OutDegreeW', 'PageRank']
    )

    parser_rankings.add_argument(
        "--edgeWeight",
        nargs="+",
        default=['RR_shrunk', 'cases_event', 'RR_shrunk x cases_event']
    )

    parser_rankings.add_argument(
        "--verbose",
        action="store_true"
    )

    args = parser.parse_args()

    if args.command == "metrics":
        compute_metrics(
            file=args.file,
            pathLoad=args.pathLoad,
            pathSave=args.pathSave,
            diseasesOut=args.diseasesOut,
            edgeWeight=args.edgeWeight,
            verbose=args.verbose,
        )

    elif args.command == "rankings":
        compute_rankings(
            file=args.file,
            pathLoad=args.pathLoad,
            pathSave=args.pathSave,
            metrics=args.metrics,
            edgeWeight=args.edgeWeight,
            verbose=args.verbose,
        )


if __name__ == "__main__":
    main()