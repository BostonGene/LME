# Python3.7
import pandas as pd
import numpy as np

from lme.classification import KNeighborsClusterClassifier
from lme.utils import read_gene_sets, ssgsea_formula, median_scale


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Classification')
    parser.add_argument('refsign', type=str,
                        help='Reference signatures')
    parser.add_argument('refannot',  type=str,
                        help='Reference annotation, Should contain MFP column')
    parser.add_argument('gmt',  type=str,
                        help='Signatures')
    parser.add_argument('exp', type=str,
                        help='Expression matrix to process')
    parser.add_argument('result',  type=str,
                        help='File to save labels')
    args = parser.parse_args()

    # Example of how to classify another cohort having clusters on TCGA for example
    # Load Training cohort with known MFP labels
    cohort_signature_scores_scaled = pd.read_csv(args.refsign, sep='\t', index_col=0).T  # Signatures in rows
    print(f'Reference signatures provided for {len(cohort_signature_scores_scaled)} samples')
    cohort_annotation = pd.read_csv(args.refannot, sep='\t', index_col=0)  # Contains MFP cluster labels in MFP column
    print(f'Reference annotation provided for {len(cohort_signature_scores_scaled)} samples')
    #  Fit the model
    MODEL = KNeighborsClusterClassifier(norm=False, scale=False, clip=2, k=35).fit(cohort_signature_scores_scaled,
                                                                                   cohort_annotation.MFP)

    # Load the cohort of interest
    # Read signatures
    gmt = read_gene_sets(args.gmt)  # GMT format like in MSIGdb
    print(f'Loaded {len(gmt)} signatures')

    # Read expressions
    exp = pd.read_csv(args.exp, sep='\t', header=0, index_col=0).T  # log2+1 transformed; Genes should appear to be in rows

    print(f'Classifying cohort, N={len(exp)} samples')
    if exp.max().max() > 35:
        print('Performing log2+1 transformation')
        exp = np.log2(1+exp)


    # Calc signature scores
    signature_scores = ssgsea_formula(exp, gmt)

    # TODO add progeny


    # Scale signatures
    signature_scores_scaled = median_scale(signature_scores, 2)

    # Predict clusters
    cluster_labels = MODEL.predict(signature_scores_scaled[MODEL.X.columns]).rename('MFP')
    print('Predicted labels count:')
    print(cluster_labels.value_counts())

    # Output the clusters
    cluster_labels.to_csv(args.result, sep='\t', index=True)
