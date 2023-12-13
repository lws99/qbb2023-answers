#!/usr/bin/env python

import sys

import scanpy as sc
import numpy
import matplotlib.pyplot as plt


#from livecoding.py import main



# Read the 10x dataset filtered down to just the highly-variable genes
adata = sc.read_h5ad("variable_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 




# EXERCISE 1 =============================================================================================================================

#Step 1.1: Computing a neighborhood graph

#Before running clustering, scanpy needs to construct what is called a “neighborhood” graph, which essentially records each cell’s “nearest neighbors”: the cells whose 
#expression looks the most similar to the focal cell.

#Using the scanpy documentation, find a pre-processing function that computes a neighborhood graph from adata. Run the function you found using the 
#parameters n_neighbors=10 and n_pcs=40.

def main():
    sc.settings.verbosity = 3
    sc.logging.print_header()#know where stuff is coming from
    adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/',
                            var_names='gene_symbols', cache=True)          #load from a matrix of data, where the data is 10x, get the file and the gene/variable names -- 2700 cells, 32738 genes

    adata.var_names_make_unique()       #make sure all genes have unique names

    sc.tl.pca(adata, svd_solver='arpack')  #doing the PCA analysis
    raw = adata.copy()
    #sc.pl.pca(adata, title='Unfiltered', save="_unfiltered.pdf")  #pl is the function for plot, all plotting needs to be under pl, gives 2 clusters, but unsure if they are meaningful, save= saves the image, show=FALSE will stop the plot from showing up again

    print("# cells, # genes before filtering:", adata.shape)
    sc.pp.filter_cells(adata, min_genes=200)            #filter out cells with really low complexity (only a few genes with expression--> most genes have 0 expression in the cells) ----- pp means preprocessed (sc.pp. is the filtering and normalizing function)
    sc.pp.filter_genes(adata, min_cells=3)              #filter out lowly informative genes
    print("# cells, # genes after filtering:", adata.shape)

    adata.var['mt'] = adata.var_names.str.startswith('MT-') #which genes are mitochondrial? --> boolean indicator of if the gene starts with MT (MT = mitochondrial gene) --> use this info to look at cell metrics (ex: what percentage of genes are mt)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)    #log1p --> do you want to do a transformation? inplace=True is very important (sc.pp.calculate_qc_metrics calculates quality control metrics)
    #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
     #        jitter=0.4, multi_panel=True)  #plot all count data using violin plots --> diff distrubutions of genes in each plot, all cells have ~2,000 transcripts --> ~2.5% of transcripts are mt
    
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]   #filter by gene counts, keeping cells with less than 2500 genes by count, get rid of cells with too much mt expression (keep cells with under 5% mt expression) --> now down to 2638 cells after removing cells with a lot of mt gene expression
    adata = adata[adata.obs.pct_counts_mt < 5, :]      
    print("# cells, # genes after MT filtering:", adata.shape)
    #sc.pl.highest_expr_genes(adata, n_top=20)

    sc.pp.normalize_total(adata, target_sum=1e4)  #normalize the expression of genes by total number of counts per cell (plot distrubution of expression across the cells want mean=0 std=1, so we can compare 2 genes and know they are simillary distributed) --> want every gene to have 10,000 transcripts expressed when normalized. To remove the 0's (0 expression), set the 0's to be 1 (all genes have at least 1 transcript) --> this allows you to log transform the expression data so you can more easily visualize it
    sc.pp.log1p(adata) #log transform expression

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
                                min_disp=0.5) #identify highly variable genes (genes that appear a lot in some cells and not others) --> min_mean expression, min_max expression, min_dispersion is measure of deviation or variance (black = highly variable, gray = not highly variable)
    #sc.pl.highly_variable_genes(adata)
    adata.write("filtered_data.h5")



    #look at PCA after everything
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    print("# cells, # genes after variability filtering:", adata.shape)

    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver='arpack')
    #sc.pl.pca_variance_ratio(adata, log=True)
    #sc.pl.pca(adata, color='CST3')

    #fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    #sc.pl.pca(raw, ax=ax[0], title="Uniltered", show=False)
    #sc.pl.pca(adata, ax=ax[1], title="Filtered", show=False)
    #plt.tight_layout()
    #plt.savefig("pca.pdf")
    #plt.close()

    adata.write('variable_data.h5')





    #EXERCISE 1=======================================


    #NEIGHBORHOOD GRAPH
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

    #Leiden clustering
    sc.tl.leiden(adata)


    #VISUALIZING CLUSTERS

    #UMAP
    sc.tl.umap(adata, maxiter=900)

    #tSNE
    sc.tl.tsne(adata)
    
    #plot clusters

    fig, axes = plt.subplots(ncols=2)

    sc.pl.umap(adata, color='leiden', ax=axes[0], title="UMAP Clusters", show=False)

    sc.pl.tsne(adata, color='leiden', ax=axes[1], title="t-SNE Clusters", show = False)



    #EXERCISE 2 =======================================================================================

    #ranking genes in each cluster

    #Wilcoxon rank-sum method
    wilcoxon_adata = sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=True, copy=True)
    print(wilcoxon_adata)


    #rank marker genes using logistic regression
    logreg_adata = sc.tl.rank_genes_groups(adata, groupby='leiden', method='logreg', use_raw=True, copy=True)
    print(logreg_adata)


    #visualize the markers


    # number in title refers to the cluster and genes with the highest score are most likely to be the marker gene for that cluster

    sc.pl.rank_genes_groups(wilcoxon_adata, n_genes=25, sharey = False, show=False, use_raw=True, save="_wilcoxon.png")


    sc.pl.rank_genes_groups(logreg_adata, n_genes=25, sharey = False, show=False, use_raw=True, save="_logreg.png")









    #EXERCISE 3 =============================================================================================

    #Reload missing genes


    leiden = adata.obs['leiden']
    umap = adata.obsm['X_umap']
    tsne = adata.obsm['X_tsne']
    adata = sc.read_h5ad('filtered_data.h5')
    adata.obs['leiden'] = leiden
    adata.obsm['X_umap'] = umap
    adata.obsm['X_tsne'] = tsne

    adata.write('filtered_clustered_data.h5')

    adata = sc.read_h5ad("filtered_clustered_data.h5")
    adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 


    #Matching genes to cell types

    #logistic genes high hits
    #LTB
    #IL32
    #CCR7
    #S100AB
    #CD79A
    #CCL5
    #NKG7
    #GZMK
    #FCGR3A
    #GNLY
    #FCER1A
    #PPBP
    

    fig, axes = plt.subplots(ncols=5)

    sc.pl.tsne(adata, ax=axes[0], title="t-SNE Clusters, PPBP Marker Gene", show = False, color = "PPBP") #best marker (marks/associated with one cluster, and strongly associated), matches expectations --only shows up for cluster 8 in logreg png
    sc.pl.tsne(adata, ax=axes[1], title="t-SNE Clusters, CCR7 Marker Gene", show = False, color = "CCR7") #ok marker (marks/associated with multiple clusters and weakly associated), matches expectations -- associaton is weak based on logreg png, even for the most associated cluster (only shows up in most associated cluster top 25 too)
    sc.pl.tsne(adata, ax=axes[2], title="t-SNE Clusters, IL32 Marker Gene", show = False, color = "IL32") # horrible marker (marks/associated with many clusters and strongly with many clusters), top 25 for 1 and 4, shows up in other clusters, but not in top 25  
    sc.pl.tsne(adata, ax=axes[3], title="t-SNE Clusters, FCER1A Marker Gene", show = False, color = "FCER1A") #another very good marker for cluster 7, matches expectations based on logreg png
    sc.pl.tsne(adata, ax=axes[4], title="t-SNE Clusters, GNLY Marker Gene", show = False, color = "GNLY") # very good marker for cluster 6
    plt.show()

main()












