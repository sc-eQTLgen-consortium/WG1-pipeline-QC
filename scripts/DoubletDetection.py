import numpy as np
import doubletdetection
import tarfile
import matplotlib.pyplot as plt
import os
import argparse
import read10x


parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix.mtx")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory; default is current working directory")
parser.add_argument("-i", "--n_iterations", required = False, default = 50, help = "Number of iterations to use; default is 50")
parser.add_argument("-p", "--phenograph", required = False, default = False, help = "Whether to use phenograph (True) or not (False); default is False")
parser.add_argument("-s", "--standard_scaling", required = False, default = True, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True")
parser.add_argument("-t", "--p_thresh", required = False, default = 1e-16, help = "P-value threshold for doublet calling; default is 1e-16")
parser.add_argument("-v", "--voter_thresh", required = False, default = 1e-16, help = "Voter threshold for doublet calling; default is 0.5")
args = parser.parse_args()

### Read in data ###
raw_counts = read10x.import_cellranger_mtx(args.counts_matrix_dir)
print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]

clf = doubletdetection.BoostClassifier(n_iters=args.n_iterations, use_phenograph=args.phenograph, standard_scaling=args.standard_scaling)
doublets = clf.fit(raw_counts).predict(p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)
np.savetxt(os.path.join(args.outdir,'DoubletDetection_predicted_doublet.txt'), doublets, fmt='%s')

### Figures ###
doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=True, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

f2, umap_coords = doubletdetection.plot.umap_plot(raw_counts, doublets, random_state=1, 
                                                       save=os.path.join(args.outdir,'umap_test.pdf'), show=True)
f3 = doubletdetection.plot.threshold(clf, save=os.path.join(args.outdir,'threshold_test.pdf'), show=True, p_step=6)