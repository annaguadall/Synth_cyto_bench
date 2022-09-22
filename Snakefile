######################################################################################
## SNAKEMAKE FILE
## SYNTHETIC SAMPLES, phenotype model_01
######################################################################################

import os
import glob

## ----Parameters --------------------------------------------------------------------
configfile: 'config.yaml'

PROJECT = config['project']

SAMPLE_LIST = config['sample_list']       # Converting a parameter of variable size as input
if not os.path.exists('sample_list'):
    os.makedirs('sample_list')
for code in SAMPLE_LIST.keys():
    if not os.path.exists('sample_list/code_{sample}.txt'.format(sample = code)):
        shell('echo {size} > sample_list/code_{sample}.txt'.format(sample = code, size = SAMPLE_LIST[code]))
'''
SAMPLES = SAMPLE_LIST.keys()
'''
SAMPLES = ['01', '02', '03', '04']
        
REPLICATES_LIST = config['replicate_list']        # Converting a parameter of variable size as input
if not os.path.exists('replicate_list'):
    os.makedirs('replicate_list')
for letter in REPLICATES_LIST.keys():
    if not os.path.exists('replicate_list/{replicate}.txt'.format(replicate = letter)):
        shell('echo {seed} > replicate_list/{replicate}.txt'.format(replicate = letter, seed = REPLICATES_LIST[letter]))
REPLICATES = REPLICATES_LIST.keys()

if config['flowsom_exact']:             
    CLUSTERS = config['flowsom_clusters']        # Converting a parameter (clusters) of variable size as input
    if not os.path.exists('flowsom_clusters'):
        os.makedirs('flowsom_clusters')
    for i in CLUSTERS:
        if not os.path.exists('flowsom_clusters/{clusters}.txt'.format(clusters = i)):
            shell('echo {clusters} > flowsom_clusters/{clusters}.txt'.format(clusters = i))

INPUT_LIST = []
if config['flowsom_max']:
    INPUT_LIST.append('results/flowsom_max/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['flowsom_exact']:
    for i in CLUSTERS:
        INPUT_LIST.append('results/flowsom_exact_' + str(i) + '/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['phenograph']:
    INPUT_LIST.append('results/phenograph/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['flowmeans']:
    INPUT_LIST.append('results/flowmeans/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['flowpeaks']:
    INPUT_LIST.append('results/flowpeaks/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['depeche']:
    INPUT_LIST.append('results/depeche/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['tsne_phenograph']:
    INPUT_LIST.append('results/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['tsne_flowmeans']:
    INPUT_LIST.append('results/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['tsne_flowpeaks']:
    INPUT_LIST.append('results/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['tsne_depeche']:
    INPUT_LIST.append('results/tsne_depeche/' + PROJECT + '_{sample}_{replicate}_summary.csv') 
if config['tsne_clusterx']:
    INPUT_LIST.append('results/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}_summary.csv') 
if config['umap_phenograph']:
    INPUT_LIST.append('results/umap_phenograph/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['umap_flowmeans']:
    INPUT_LIST.append('results/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['umap_flowpeaks']:
    INPUT_LIST.append('results/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}_summary.csv')
if config['umap_depeche']:
    INPUT_LIST.append('results/umap_depeche/' + PROJECT + '_{sample}_{replicate}_summary.csv') 
if config['umap_clusterx']:
    INPUT_LIST.append('results/umap_clusterx/' + PROJECT + '_{sample}_{replicate}_summary.csv') 


## ----Rule all ----------------------------------------------------------------------
rule all:
    input:
        expand(INPUT_LIST, sample = SAMPLES, replicate = REPLICATES),
        'results/summary.txt',
        expand('exprs_plots/' + PROJECT + '_{sample}_{replicate}_density_CD4.png', sample = SAMPLES, replicate = REPLICATES)

## ----Generation of synthetic samples -----------------------------------------------        
rule generate_fcs:
    input:
        sample = 'sample_list/code_{sample}.txt',
        replicate = 'replicate_list/{replicate}.txt'
    params:
        neg_mean = config['neg_mean'],
        neg_sd = config['neg_sd'],
        pos_mean = config['pos_mean'],
        pos_sd = config['pos_sd'],
        phenotypes = config['phenotypes']
    output:
        fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs',
        labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        
    log:
        'logs/generate/' + PROJECT + '_{sample}_{replicate}.log'
    script:
        'scripts/pop_fcs_generator.R'
    
## ----Plotting fcs ----------------------------------------------------------------       
include: config['plot_fcs_snakefile']


## ----Reduction to 2 dimensions ---------------------------------------------------
if config['tsne']:
    rule tsne:
        input:
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs'
        output:
            red = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne.fcs',
            runtime = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne_time.csv'
        log:
            'logs/tsne/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/tsne.R'
 
if config['umap']:
    rule umap:
        input:
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs'
        output:
            red = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap.fcs',
            runtime = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap_time.csv'
        log:
            'logs/umap/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/umap.R'   
    
            
## ----Benchmarking ----------------------------------------------------------------

## Clustering on all dimensions

if config['flowsom_max']:
    rule flowsom_max:
        input:
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/flowsom_max/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/flowsom_max/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/flowsom_max/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/flowsom_max/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/flowsom_max/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/flowsom_max/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'flowsom_max',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/flowsom_max/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowsom_max.R'

if config['flowsom_exact']:             
    rule flowsom_exact:
        input:
            n_clus = 'flowsom_clusters/{clusters}.txt',
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv',
        output:
            clus = 'results/flowsom_exact_{clusters}/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/flowsom_exact_{clusters}/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/flowsom_exact_{clusters}/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/flowsom_exact_{clusters}/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/flowsom_exact_{clusters}/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/flowsom_exact_{clusters}/' + PROJECT + '_{sample}_{replicate}_summary.csv',
        params:
            method = 'flowsom_exact_{clusters}',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/flowsom_exact_{clusters}/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowsom_exact.R'

if config['phenograph']:
    rule phenograph:
        input:
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/phenograph/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/phenograph/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/phenograph/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/phenograph/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/phenograph/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/phenograph/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'phenograph',
            reduction_time = 'FALSE',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/phenograph/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/phenograph.R'

if config['flowmeans']:
    rule flowmeans:
        input:
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/flowmeans/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/flowmeans/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/flowmeans/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/flowmeans/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/flowmeans/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/flowmeans/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'flowmeans',
            reduction_time = 'FALSE',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/flowmeans/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowmeans.R'

if config['flowpeaks']:
    rule flowpeaks:
        input:
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/flowpeaks/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/flowpeaks/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/flowpeaks/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/flowpeaks/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/flowpeaks/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/flowpeaks/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'flowpeaks',
            reduction_time = 'FALSE',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/flowpeaks/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowpeaks.R'

if config['depeche']:
    rule depeche:
        input:
            fcs = 'data/' + PROJECT + '_{sample}_{replicate}.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/depeche/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/depeche/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/depeche/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/depeche/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/depeche/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/depeche/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'depeche',
            reduction_time = 'FALSE',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/depeche/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/depeche.R'
            

## Clustering on t-SNE-reduced data

if config['tsne_phenograph']:
    rule tsne_phenograph:
        input:
            fcs = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'tsne_phenograph',
            reduction_time = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/tsne_phenograph/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/phenograph.R'

if config['tsne_flowmeans']:
    rule tsne_flowmeans:
        input:
            fcs = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'tsne_flowmeans',
            reduction_time = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/tsne_flowmeans/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowmeans.R'

if config['tsne_flowpeaks']:
    rule tsne_flowpeaks:
        input:
            fcs = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'tsne_flowpeaks',
            reduction_time = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/tsne_flowpeaks/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowpeaks.R'

if config['tsne_depeche']:
    rule tsne_depeche:
        input:
            fcs = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/tsne_depeche/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/tsne_depeche/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/tsne_depeche/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/tsne_depeche/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/tsne_depeche/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/tsne_depeche/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'tsne_depeche',
            reduction_time = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/tsne_depeche/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/depeche.R'
            
if config['tsne_clusterx']:
    rule tsne_clusterx:
        input:
            fcs = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'tsne_clusterx',
            reduction_time = 'reduced/tsne/' + PROJECT + '_{sample}_{replicate}_tsne_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/tsne_clusterx/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/clusterx.R'            
            
            
## Clustering on UMAP-reduced data           
            
if config['umap_phenograph']:
    rule umap_phenograph:
        input:
            fcs = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/umap_phenograph/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/umap_phenograph/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/umap_phenograph/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/umap_phenograph/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/umap_phenograph/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/umap_phenograph/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'umap_phenograph',
            reduction_time = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/umap_phenograph/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/phenograph.R'

if config['umap_flowmeans']:
    rule umap_flowmeans:
        input:
            fcs = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'umap_flowmeans',
            reduction_time = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/umap_flowmeans/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowmeans.R'

if config['umap_flowpeaks']:
    rule umap_flowpeaks:
        input:
            fcs = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'umap_flowpeaks',
            reduction_time = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/umap_flowpeaks/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/flowpeaks.R'

if config['umap_depeche']:
    rule umap_depeche:
        input:
            fcs = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/umap_depeche/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/umap_depeche/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/umap_depeche/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/umap_depeche/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/umap_depeche/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/umap_depeche/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'umap_depeche',
            reduction_time = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/umap_depeche/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/depeche.R'
            
if config['umap_clusterx']:
    rule umap_clusterx:
        input:
            fcs = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap.fcs',
            labels = 'data/' + PROJECT + '_{sample}_{replicate}_labels.csv'
        output:
            clus = 'results/umap_clusterx/' + PROJECT + '_{sample}_{replicate}_clustering.csv',
            cont = 'results/umap_clusterx/' + PROJECT + '_{sample}_{replicate}_contingency_matrix.csv',
            cm = 'results/umap_clusterx/' + PROJECT + '_{sample}_{replicate}_confusion_matrix.csv',
            cmm = 'results/umap_clusterx/' + PROJECT + '_{sample}_{replicate}_confusion_matrix_merged.csv',
            f1m = 'results/umap_clusterx/' + PROJECT + '_{sample}_{replicate}_f1_matrix.csv',
            summ = 'results/umap_clusterx/' + PROJECT + '_{sample}_{replicate}_summary.csv'
        params:
            method = 'umap_clusterx',
            reduction_time = 'reduced/umap/' + PROJECT + '_{sample}_{replicate}_umap_time.csv',
            synth = config['synthetic'],
            thres = config['threshold'],
            norm = config['normalization'],
            sample = PROJECT + '_{sample}_{replicate}'
        log:
            'logs/umap_clusterx/' + PROJECT + '_{sample}_{replicate}.log'
        script:
            'scripts/clusterx.R'          


## ----Summary ---------------------------------------------------------------------

# Summary
rule summary:
    input:
        expand(INPUT_LIST, sample = SAMPLES, replicate = REPLICATES)
    output:
        'results/summary.txt'
    shell:
        "for file in {input}; do cat $file | tail -n 1 | cut -d ',' -f2- >> {output}; done"

# Adding a header
HEADER = 'File, Norm, n, Method, Threshold, Clusters, Partitions, Mean_F1, Weighted_mean_F1, \
         Inversed_weights_mean_F1, Mean_F1_merged, Corrected_mean_F1_merged, Reduction_user_time, \
         Reduction_elapsed_time, Clustering_user_time, Clustering_elapsed_time'
         
shell('echo {header} > results/{project}_summary.txt'.format(header = HEADER, project = PROJECT))
shell('cat results/summary.txt >> results/{project}_summary.txt'.format(header = HEADER, project = PROJECT))
