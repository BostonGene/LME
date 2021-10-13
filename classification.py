# Python3.7
import pandas as pd
import numpy as np

from lme.classification import KNeighborsClusterClassifier
from lme.pathway_scoring import run_progeny
from lme.utils import read_gene_sets, ssgsea_formula, median_scale, to_common_samples

signatures_selected = [
    'Lymphatic_endothelium',
    'Angiogenesis',
    'CAF',
    'Fibroblastic_reticular_cells',
    'Matrix',
    'Matrix_remodeling',
    'Granulocyte_traffic',
    'Protumor_cytokines',
    'Follicular_dendritic_cells',
    'Macrophages',
    'M1_signature',
    'T_cell_traffic',
    'MHCII',
    'MHCI',
    'Follicular_B_helper_T_cells',
    'Treg',
    'T_cells',
    'Checkpoint_inhibition',
    'NK_cells',
    'B_cells_traffic',
    'B_cells',
    'Proliferation_rate']

progeny_selected = ['NFkB', 'p53', 'PI3K']


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

    cohort_ann_filtered = cohort_annotation[(cohort_annotation.Diagnosis == 'Diffuse_Large_B_Cell_Lymphoma') &
                                            (~cohort_annotation.LME.isna())]
    print(f'Using {len(cohort_ann_filtered)} DLBCL samples with known LME status')

    LME_MODEL = KNeighborsClusterClassifier(norm=False, scale=False, clip=3, k=35).fit(
        *to_common_samples([cohort_signature_scores_scaled[signatures_selected + progeny_selected],
                            cohort_ann_filtered.LME]))

    # Load the cohort of interest
    # Read signatures
    gmt = read_gene_sets(args.gmt)  # GMT format like in MSIGdb
    print(f'Loaded {len(gmt)} signatures, using {len(signatures_selected)} selected')

    # Read expressions
    exp = pd.read_csv(args.exp, sep='\t', header=0, index_col=0).T  # log2+1 transformed; Genes should appear to be in rows

    print(f'Classifying cohort, N={len(exp)} samples')
    if exp.max().max() > 35:
        print('Performing log2+1 transformation')
        exp = np.log2(1+exp)

    # Calc signature scores both ssgsea and progeny
    signature_scores = pd.concat([ssgsea_formula(exp, gmt), run_progeny(exp)], axis=1)

    # Scale signatures
    signature_scores_scaled = median_scale(signature_scores, 2)

    # Predict clusters
    cluster_labels = LME_MODEL.predict(signature_scores_scaled[LME_MODEL.X.columns]).rename('LME')
    print('Predicted labels count:')
    print(cluster_labels.value_counts())

    # Output the clusters
    cluster_labels.to_csv(args.result, sep='\t', index=True)
