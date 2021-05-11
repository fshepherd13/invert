import pandas as pd

#Ratios of mRNA:total mRNA for PB2, PB1, PA, HA, NP, NA, M AND S
ratios = pd.read_table("../invert_results/4_12hr_Ciliated_4_S4_sub/cmrna/4_12hr_Ciliated_4_S4_sub_cmratio.txt", index_col=0, keep_default_na=False)

splice_counts = pd.read_table("../invert_results/4_12hr_Ciliated_4_S4_sub/splice_counts/4_12hr_Ciliated_4_S4_sub_splicing_count.txt")

expression_levels = pd.read_table("../cufflinks/4_12hr_Ciliated_4_S4_sub/genes.fpkm_tracking.txt", index_col=4, keep_default_na=False)

############First, perform calculations for M gene############
#Define c(M2) <-read counts of spliced M2
cM2 = (1/2)*((splice_counts['Depth_total_left'].values[0] - splice_counts['Depth_unspliced_left'].values[0]) + (splice_counts['Depth_total_right'].values[0] - splice_counts['Depth_unspliced_right'].values[0]))

#Define c(mM) <- read counts of total mRNA of M gene
mrna_total_rna_ratio_m = ratios.loc['M', 'mRNA:total_pos_RNA']

cmM = (1/2)*(splice_counts['Depth_total_left'].values[0] + splice_counts['Depth_total_right'].values[0]) * mrna_total_rna_ratio_m

#define "f_M" <- fraction of spliced M2 read counts to total M mRNA read counts:
f_M = cM2 / cmM

############Perform same calculations for NS############
cNEP = (1/2)*((splice_counts['Depth_total_left'].values[1] - splice_counts['Depth_unspliced_left'].values[1]) + (splice_counts['Depth_total_right'].values[1] - splice_counts['Depth_unspliced_right'].values[1]))

mrna_total_rna_ratio_ns = ratios.loc['NS', 'mRNA:total_pos_RNA']
cmNS = (1/2)*(splice_counts['Depth_total_left'].values[1] + splice_counts['Depth_total_right'].values[1]) * mrna_total_rna_ratio_ns

f_NS = cNEP / cmNS

#Pull FPKM values for positive sense transcripts
HA_fpkm = expression_levels.loc['HA', 'FPKM']
NA_fpkm = expression_levels.loc['NA', 'FPKM']
M_fpkm = expression_levels.loc['cM_mM1,mM2', 'FPKM']
NP_fpkm = expression_levels.loc['NP', 'FPKM']
PA_fpkm = expression_levels.loc['PA', 'FPKM']
NS_fpkm = expression_levels.loc["cNS_mNS1,mNEP", 'FPKM']
PB1_fpkm = expression_levels.loc["PB1", 'FPKM']
PB2_fpkm = expression_levels.loc["PB2", 'FPKM']

#Pull FPKM values for negative sense transcripts
vHA_fpkm = expression_levels.loc['vHA', 'FPKM']
vNA_fpkm = expression_levels.loc['vNA', 'FPKM']
vM_fpkm = expression_levels.loc['vM', 'FPKM']
vNP_fpkm = expression_levels.loc['vNP', 'FPKM']
vPA_fpkm = expression_levels.loc['vPA', 'FPKM']
vNS_fpkm = expression_levels.loc['vNS', 'FPKM']
vPB1_fpkm = expression_levels.loc['vPB1', 'FPKM']
vPB2_fpkm = expression_levels.loc['vPB2', 'FPKM']

#Create dataframe to hold all results
#First, insert the FPKM values from the Cufflinks results and mrna:(mrna+crna) ratios
results = pd.DataFrame(
    [[HA_fpkm, vHA_fpkm, ratios.loc['HA', 'mRNA:total_pos_RNA']],
    [NA_fpkm, vNA_fpkm, ratios.loc['NA', 'mRNA:total_pos_RNA']],
    [M_fpkm, vM_fpkm, ratios.loc['M', 'mRNA:total_pos_RNA']],
    [NP_fpkm, vNP_fpkm, ratios.loc['NP', 'mRNA:total_pos_RNA']],
    [PA_fpkm, vPA_fpkm, ratios.loc['PA', 'mRNA:total_pos_RNA']],
    [NS_fpkm, vNS_fpkm, ratios.loc['NS', 'mRNA:total_pos_RNA']],
    [PB1_fpkm, vPB1_fpkm, ratios.loc['PB1', 'mRNA:total_pos_RNA']],
    [PB2_fpkm, vPB2_fpkm, ratios.loc['PB2', 'mRNA:total_pos_RNA']],
    [None, None, None], 
    [None, None, None],
    [None, None, None],
    [None, None, None]],
    columns=['FPKM_positive_sense', 'FPKM_negative_sense', 'mrna:(mrna,crna)_ratio'],
    index = ['HA', 'NA', 'M', 'NP', 'PA', 'NS', 'PB1', 'PB2', 'NS1', 'NEP', 'M1', 'M2'])

total_FPKM = sum([HA_fpkm, NA_fpkm, M_fpkm, NP_fpkm, PA_fpkm, NS_fpkm, PB1_fpkm, PB2_fpkm, vHA_fpkm, vNA_fpkm, vM_fpkm, vNP_fpkm, vPA_fpkm, vNS_fpkm, vPB1_fpkm, vPB2_fpkm])

#Transform FPKM to TPM for positive and negative sense RNA transcripts
results['TPM_positive_sense'] = results['FPKM_positive_sense'] / total_FPKM
results['TPM_negative_sense'] = results['FPKM_negative_sense'] / total_FPKM

#Create column for vRNA TPM counts (same as the direct cufflinks output)
results['vrna_TPM'] = results['TPM_negative_sense']

#Calculate crna TPM
results['crna_TPM'] = results['TPM_positive_sense'] / (1-results['mrna:(mrna,crna)_ratio'])


#mrna = pd.DataFrame([results.loc['HA', 'TPM_positive_sense'] * results.loc['HA', 'mrna:(mrna,crna)_ratio']]])

def mrna_calc(gene):
    mrna = results.loc[str(gene), 'TPM_positive_sense'] * results.loc[str(gene), 'mrna:(mrna,crna)_ratio']
    return mrna

#initialize empty dictionary to hold majority of mrna calculations:
d = dict()
for x in ['HA', 'NA', 'NP', 'PA', 'PB1', 'PB2']:
    d[x] = mrna_calc(x)

#M and NS don't exist as mrna, so add these entries to dictionary with "NA"
d['M'] = None
d['NS'] = None


#Calculate mrna for unspliced M1 and NS1:
d['NS1'] = results.loc['NS', 'TPM_positive_sense'] * (1 - f_NS) * results.loc['NS', 'mrna:(mrna,crna)_ratio']

d['M1'] = results.loc['M', 'TPM_positive_sense'] * (1-f_M) * results.loc['M', 'mrna:(mrna,crna)_ratio']

#Calculate mrna for spliced M2 and NEP:
d['NEP'] = results.loc['NS', 'TPM_positive_sense'] * f_NS * results.loc['NS', 'mrna:(mrna,crna)_ratio']

d['M2'] = results.loc['M', 'TPM_positive_sense'] * (f_M) * results.loc['M', 'mrna:(mrna,crna)_ratio']

#Map the full dictionary back to the results dataframe
results['mrna_TPM'] = results.index.map(d)

results.to_csv("../invert_results/4_12hr_Ciliated_4_S4_sub/results.csv")