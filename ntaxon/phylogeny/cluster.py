import pandas as pd
from scipy.spatial import distance
from scipy import cluster
from math import isnan
from ntaxon.fingerprinting import BinaryMatrix


# Find out best Linkage Method having highest cophenetic correlation coefficient
def linkage_info(pdist_val, linkage_method: str) -> dict:
    # ['Linkage Method', 'Linkage Validity', 'Inconsistency', 'Monotony', 'Cophenetic Correlation']
    _linkage = cluster.hierarchy.linkage(pdist_val, method=linkage_method)
    corr, coph_dist = cluster.hierarchy.cophenet(_linkage, pdist_val)

    return {
      "Linkage": linkage_method,
      "Linkage Validity": cluster.hierarchy.is_valid_linkage(_linkage),
      "Inconsistency": cluster.hierarchy.is_valid_im(_linkage),
      "Monotony": cluster.hierarchy.is_monotonic(_linkage),
      "Copenetic Correlation": corr
    }


def screen_methods(data: BinaryMatrix) -> pd.DataFrame:
    screening_results = pd.DataFrame(columns=["Distance", "Linkage", "Linkage Validity", "Inconsistency", "Monotony", "Copenetic Correlation"])

    distance_methods = ['dice', 'hamming', 'jaccard', 'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',
                    'sokalsneath', 'yule']
    linkage_methods = ['average', 'median', 'single', 'complete', 'centroid', 'ward']

    # TODO Compute custom distance methods such as Nei and Li

    for _dm in distance_methods:
        mat = data
        if isinstance(data, BinaryMatrix):
            mat = data.get_binary_matrix()
        _pdist_val = 1 - distance.pdist(mat, _dm)

        if not any(_pdist_val < 0):  # check if any negative value in distance matrix
            for method in linkage_methods:
                _info = linkage_info(_pdist_val, method)
                if not isnan(_info['Copenetic Correlation']):
                    _info['Distance'] = _dm
                    screening_results = screening_results.append(_info, ignore_index=True)
    return screening_results.sort_values(['Copenetic Correlation'], ascending=False).reset_index(drop=True)


def optimal_clustering_params(data: BinaryMatrix, monotonic=True):
    # TODO: Fix broken
    res = screen_methods(data)
    res = res[res['Monotony'] == monotonic]
    print(res)
    print(res.iloc[1]['Distance', 'Linkage'])
    #return tuple(res.iloc[1]['Distance', 'Linkage'].to_list())
