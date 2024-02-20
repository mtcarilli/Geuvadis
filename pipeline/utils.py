# functions for normalizing gene expression and creating .bed file from gtf and gene expression matrix
import pandas as pd
import numpy as np
import scipy.stats as stats
import gzip
import functools



def build_gene_location_df(gtf_file_path, genes_to_use):
    genes = []
    chromosomes = []
    s1_list = []
    s2_list = []
    
    if gtf_file_path[-3:] == '.gz':
        open_fn = functools.partial(gzip.open, mode='rt')
    else:
        open_fn = open
    with open_fn(gtf_file_path) as gtf_file:
        for line in gtf_file:
            if line[0] == '#':  # skip comments
                continue
            l = line.strip().split()
            feat_type = l[2]
            gene_name = l[9][1:-2]
            if (feat_type) == 'gene' and (gene_name in genes_to_use):
                strand = l[6]
                chrom = str(l[0])
                if strand == '+' and chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '10', '12','13', '14', '15', '16', '17', '18', '20', '19', '21', '22']:
                    s1 = np.int64(l[3])-1
                    s2 = np.int64(l[3])
                    genes.append(gene_name)
                    chromosomes.append(str(chrom))
                    s1_list.append(s1)
                    s2_list.append(s2)
                elif strand == '-' and chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '10', '12','13', '14', '15', '16', '17', '18', '20', '19', '21', '22']:
                    s1 = np.int64(l[4])-1
                    s2 = np.int64(l[4])
                    genes.append(gene_name)
                    chromosomes.append(str(chrom))
                    s1_list.append(s1)
                    s2_list.append(s2)
                # elif chrom not in ['X', 'Y', 'M']:
                #     s1 = l[3]
                #     s2 = l[4]
                #     genes.append(gene_name)
                #     chromosomes.append(chrom)
                #     s1_list.append(s1)
                #     s2_list.append(s2)
                    
    return pd.DataFrame({'chr' : chromosomes, 'start' : s1_list, 'end' : s2_list,'gene_id' : genes})


def get_gene_lengths_df(gtf_file_path, genes_to_use):
    genes = []
    chromosomes = []
    s1_list = []
    s2_list = []
    
    if gtf_file_path[-3:] == '.gz':
        open_fn = functools.partial(gzip.open, mode='rt')
    else:
        open_fn = open
    with open_fn(gtf_file_path) as gtf_file:
        for line in gtf_file:
            if line[0] == '#':  # skip comments
                continue
            l = line.strip().split()
            feat_type = l[2]
            gene_name = l[9][1:-2]
            if (feat_type) == 'gene' and (gene_name in genes_to_use):
                strand = l[6]
                chrom = str(l[0])
                if strand == '+' and chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '10', '12','13', '14', '15', '16', '17', '18', '20', '19', '21', '22']:
                    s1 = np.int64(l[3])
                    s2 = np.int64(l[4])
                    genes.append(gene_name)
                    chromosomes.append(str(chrom))
                    s1_list.append(s1)
                    s2_list.append(s2)
                elif strand == '-' and chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '10', '12','13', '14', '15', '16', '17', '18', '20', '19', '21', '22']:
                    s1 = np.int64(l[3])
                    s2 = np.int64(l[4])
                    genes.append(gene_name)
                    chromosomes.append(str(chrom))
                    s1_list.append(s1)
                    s2_list.append(s2)
                # elif chrom not in ['X', 'Y', 'M']:
                #     s1 = l[3]
                #     s2 = l[4]
                #     genes.append(gene_name)
                #     chromosomes.append(chrom)
                #     s1_list.append(s1)
                #     s2_list.append(s2)

    gene_order_df_ = pd.DataFrame({'gene_id' : genes_to_use})
    gene_length_df_ = pd.DataFrame({'chr' : chromosomes, 'start' : s1_list, 'end' : s2_list,'gene_id' : genes})
    gene_lengths_df = gene_order_df_ .merge(gene_lengths_df_, on = 'gene_id')
    return gene_lengths_df

def get_phenotype_bed(adata,layer,gtf_path,save_path):
  '''saves phenotype values in .bed format to save path.

  Format:
  #chr start end gene_id samp1 .... sampN
  1    100   101 ENSGXXXX 0.0 ..... 45.8

  adata: adata objectsubset to desired genes
  layer: layer of the adata object to use as values
  gtf_path: path to gtf to use
  '''

  genes_to_use = [g.split('.')[0] for g in adata.var.var_names.values]
  gene_df = build_gene_location_df(gtf_path,genes_to_use=genes_to_use)

  adata.var['gene_id'] = genes_to_use
  adata.obs.index = adata.obs.obs_names
  expression_df = pd.DataFrame( {ind : np.array(adata[adata.obs.Individual==ind].layers[layer].todense()).flatten() for ind in adata.obs.Individual.values} ) 
  expression_df['gene_id'] = adata.var.gene_id.values
  bed_df = pd.merge(gene_df,expression_df,on='gene_id',sort=False)
  bed_df = bed_df.sort_values(['chr', 'start'], ascending=[True, True])
  bed_df.to_csv(save_path,sep='\t',index=False)
  # !sed -i '1s/^/#/' $save_path

  
# functions for normalization 
def quantile_normalize(df):
    """Quantile normalizes a pandas DataFrame.
      From PCCA github
      https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe
      Args:
        df: A pandas DataFrame
      Returns:
        The dataframe where all columns have the same distribution.
      """
    rank_mean = df.stack().groupby(
    df.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

# useful function for unpacking arrays
def nd(arr):
    return np.asarray(arr).reshape(-1)


def inverse_normal_transform(array):
    ''' Transform each row of expression (sample expression) to the inverse normal distribution. 
    USES BLOM METHOD: https://www.statsdirect.com/help/data_preparation/transform_normal_scores.htm 
    '''

    array_T = array.T # want to transform the row, being genes, to inverse normal
    
    # rank each array 
    new_array = np.ones_like(array,dtype = float)
    for i,row in enumerate(array):
        for j,val in enumerate(row):
            val_quantile =  (sum(row<val) + .375)/(len(row) + 1.4)
            inv_val = stats.norm.ppf(val_quantile)
            new_array[i,j] = inv_val
            
    return(new_array.T)

def TPM(array,gene_lengths):
    ''' 
    1. Divide read counts by length of gene in kilobases --> RPK
    2. Sum all RPK and divide by 1,000,000 --> per million scaling factor for given sample.
    3. Divide RPK by scaling factor --> TPM.
    '''
    gene_lengths = gene_lengths/1000.0
    rpk_array = array/gene_lengths
    scaling_factor = rpk_array.sum(axis=1)/1000000.0
    tpm_array = rpk_array/scaling_factor[:,None]

    return(tpm_array)

def RPKM(array,gene_lengths):
    ''' 
    1. Sum all RPK and divide by 1,000,000 --> per million scaling factor for given sample.
    2. Divide all reads by scaling factor --> reads per million RPM.
    3. Divide RPM by length of gene in kilobases --> RPKM.
    '''

    gene_lengths = gene_lengths/1000.0
    scaling_factor = array.sum(axis=1)/1000000.0
    rpm_array = array/scaling_factor[:,None]
    rpkm_array = rpm_array/gene_lengths

    return(rpkm_array)

def log1p(array):
     return(np.log1p(array))





