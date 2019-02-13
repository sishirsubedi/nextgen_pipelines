import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
warnings.filterwarnings("ignore")

SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
ENV=sys.argv[3]

### filter parameters
VARSCAN_PVAL=0.01
NORMAL_DEPTH=25
TUMOR_DEPTH=20
NORMAL_ALT_ALLELE_FREQ=2.0
TUMOR_ALT_ALLELE_FREQ=10.0
NORMAL_STRAND_BIAS_PLUS=15
NORMAL_STRAND_BIAS_MINUS=4
TUMOR_STRAND_BIAS_1_MINUS=0
TUMOR_STRAND_BIAS_1_PLUS=0
TUMOR_STRAND_BIAS_2_MINUS=2
TUMOR_STRAND_BIAS_2_PLUS=8

DETAIL=1


def checkRef(df_snp_mut_2,detail):
    df_snp_1 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/reference_standard_creDesign.txt",sep='\t')
    df_snp_mut_1 = df_snp_1[['CHROM','POS','REF','ALT']]

    old_cols = list(df_snp_mut_2.columns)
    old_cols[1:5] = ['CHROM','POS','REF','ALT']
    df_snp_mut_2.columns = old_cols

    df_join = pd.merge(left=df_snp_mut_1, right=df_snp_mut_2, on=['CHROM','POS','REF','ALT'],how='outer',indicator=True)
    # df_join.to_csv("temp.txt",sep='\t')

    if detail:
        sns.distplot(df_join[df_join['_merge']=='both']['normal_var_freq_percent'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['normal_var_freq_percent'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_normal_var_freq_percent')
        plt.close()

        plt.xlim(0,2.5)
        sns.distplot(df_join[df_join['_merge']=='both']['normal_var_freq_percent'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['normal_var_freq_percent'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_normal_var_freq_percent_2')
        plt.close()


        sns.distplot(df_join[df_join['_merge']=='both']['tumor_var_freq_percent'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['tumor_var_freq_percent'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_var_freq_percent')
        plt.close()

        sns.distplot(df_join[df_join['_merge']=='both']['normal_depth'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['normal_depth'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_normal_depth')
        plt.close()

        sns.distplot(df_join[df_join['_merge']=='both']['tumor_depth'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['tumor_depth'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_depth')
        plt.close()

        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['normal_depth'] < 50 ) )]['normal_depth'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['normal_depth'] < 50) )]['normal_depth'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_normal_depth_2')
        plt.close()

        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['tumor_depth'] < 50 ) )]['tumor_depth'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['tumor_depth'] < 50 ))]['tumor_depth'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_depth_2')
        plt.close()

        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['normal_reads1_plus'] < 50 ) )]['normal_reads1_plus'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['normal_reads1_plus'] < 50) )]['normal_reads1_plus'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_normal_reads1_plus')
        plt.close()


        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['normal_reads1_minus'] < 50 ) )]['normal_reads1_minus'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['normal_reads1_minus'] < 50) )]['normal_reads1_minus'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_normal_reads1_minus')
        plt.close()

        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['tumor_reads1_plus'] < 50 ) )]['tumor_reads1_plus'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['tumor_reads1_plus'] < 50) )]['tumor_reads1_plus'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads1_plus')
        plt.close()


        sns.distplot(df_join[df_join['_merge']=='both' ]['tumor_reads1_plus'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['tumor_reads1_plus'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads1_plus_2')
        plt.close()


        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['tumor_reads1_minus'] < 50 ) )]['tumor_reads1_minus'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['tumor_reads1_minus'] < 50) )]['tumor_reads1_minus'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads1_minus')
        plt.close()


        sns.distplot(df_join[df_join['_merge']=='both']['tumor_reads1_minus'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['tumor_reads1_minus'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads1_minus_2')
        plt.close()


        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['tumor_reads2_plus'] < 50 ) )]['tumor_reads2_plus'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['tumor_reads2_plus'] < 50) )]['tumor_reads2_plus'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads2_plus')
        plt.close()

        sns.distplot(df_join[df_join['_merge']=='both']['tumor_reads2_plus'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['tumor_reads2_plus'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads2_plus_2')
        plt.close()


        plt.xlim(0,50)
        sns.distplot(df_join[( (df_join['_merge']=='both') & (df_join['tumor_reads2_minus'] < 50 ) )]['tumor_reads2_minus'], color='green',kde=False,label='TP',bins=50)
        sns.distplot(df_join[(  (df_join['_merge']=='right_only') & (df_join['tumor_reads2_minus'] < 50) )]['tumor_reads2_minus'], color='red',kde=False,label='FP',bins=50)
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads2_minus')
        plt.close()


        sns.distplot(df_join[df_join['_merge']=='both']['tumor_reads2_minus'], color='green',kde=False,label='TP')
        sns.distplot(df_join[df_join['_merge']=='right_only']['tumor_reads2_minus'], color='red',kde=False,label='FP')
        plt.legend()
        plt.savefig(OUT_DIR+"/"+SAMPLE+'_DIST_PLOT_tumor_reads2_minus_2')
        plt.close()


    return (df_join[df_join['_merge']=='both'].shape[0],df_join[df_join['_merge']=='right_only'].shape[0])


# for FDR in [0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001,0.000000001,0.0000000001,0.00000000001]:
for FDR in [0.0000000001]:
    # print(SAMPLE+OUT_DIR+ENV)

    file_snp=OUT_DIR+SAMPLE+".varscan.output.snp"
    file_indel=OUT_DIR+SAMPLE+".varscan.output.indel"
    df_snp = pd.read_csv(file_snp,sep='\t')
    df_snp['chrom'] =df_snp['chrom'].apply(str)
    df_snp['group']='snp'

    df_indel = pd.read_csv(file_indel,sep='\t')
    df_indel['chrom'] =df_indel['chrom'].apply(str)
    df_indel['group']='indel'
    for indx, row in df_indel.iterrows():
        if row['var'][0] == '-':
            df_indel.set_value(indx, "ref", row['ref']+row['var'][1:])
            df_indel.set_value(indx, "var", row['ref'])
        elif row['var'][0] == '+':
            df_indel.set_value(indx, "var", row['ref'] + row['var'][1:] )

    df_both = pd.concat([df_snp,df_indel])
    df_both.sort_values(by=['somatic_status'],inplace=True)
    # print('Total Raw:'+ str(df_both.shape[0]))
    ## plot raw counts
    # sns.countplot(x=df_both['somatic_status'],hue=df_both['group'],data=pd.melt(df_both))
    # plt.savefig(OUT_DIR+"/"+'snp_indel_count_raw')
    # plt.close()
    # print('plot completed')

    ###filtering:
    filter_indx = df_both[df_both['chrom'].str.contains("GL|gl")].index
    df_both.drop(filter_indx, inplace=True)
    print('After GL/gl:'+ str(df_both.shape[0]))
    print('Raw:snp:'+ str(df_both[df_both['group']=='snp'].shape[0]))
    print('Raw:indel:'+ str(df_both[df_both['group']=='indel'].shape[0]))


    # somatic only
    df_both_pass_s0 = df_both[(df_both['somatic_status'].isin(['Somatic']) )]
    print('Somatic:'+ str(df_both_pass_s0.shape[0]))
    print('Somatic:snp:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='snp'].shape[0]))
    print('Somatc:indel:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='indel'].shape[0]))

    df_both_pass_s0['normal_var_freq_percent'] = df_both_pass_s0['normal_var_freq'].str.replace('%','' ).astype(float)
    df_both_pass_s0['tumor_var_freq_percent'] = df_both_pass_s0['tumor_var_freq'].str.replace('%','' ).astype(float)
    df_both_pass_s0['normal_depth'] = df_both_pass_s0['normal_reads1'] + df_both_pass_s0['normal_reads2']
    df_both_pass_s0['tumor_depth'] = df_both_pass_s0['tumor_reads1'] + df_both_pass_s0['tumor_reads2']

    def pvalAdj(p_val_rank,total,fdr):
        return (p_val_rank*fdr)/float(total)

    df_both_pass_s0.sort_values(by=['somatic_p_value'],inplace=True)
    df_both_pass_s0.reset_index(inplace=True)
    df_both_pass_s0['somatic_p_value_rank']= [ x for x in range(1,df_both_pass_s0.shape[0]+1)]
    df_both_pass_s0['somatic_p_value_adj']= [ pvalAdj(x,df_both_pass_s0.shape[0],FDR) for x in df_both_pass_s0['somatic_p_value_rank']]
    df_both_pass_s0['somatic_p_value_adj_filter']=[1 if x[0]<x[1] else 0 for x in zip(df_both_pass_s0['somatic_p_value'], df_both_pass_s0['somatic_p_value_adj']) ]
    df_both_pass_s0 = df_both_pass_s0[df_both_pass_s0['somatic_p_value_adj_filter']==1]
    print('After pvalue:'+ str(df_both_pass_s0.shape[0]))

    vals=checkRef(df_both_pass_s0,DETAIL)
    print(str(FDR)+'\t'+str(vals[0])+'\t'+str(vals[1]))

    # df_both_pass_s0.to_csv("_temp.txt",sep='\t',index=False)



    # normal var freq
    df_both_pass_s0 = df_both_pass_s0[df_both_pass_s0['normal_var_freq_percent']<NORMAL_ALT_ALLELE_FREQ]
    df_both_pass_s0 = df_both_pass_s0[df_both_pass_s0['tumor_var_freq_percent']>=TUMOR_ALT_ALLELE_FREQ]
    print('After normal var freq:'+ str(df_both_pass_s0.shape[0]))

    vals=checkRef(df_both_pass_s0,0)
    print(str(FDR)+'\t'+str(vals[0])+'\t'+str(vals[1]))


    # # min depth 20x for both normal and tumor_reads
    df_both_pass_s0 = df_both_pass_s0[ (df_both_pass_s0['normal_depth']>=NORMAL_DEPTH) & (df_both_pass_s0['tumor_depth']>=TUMOR_DEPTH) ]
    print('After depth 20x:'+ str(df_both_pass_s0.shape[0]))

    vals=checkRef(df_both_pass_s0,0)
    print(str(FDR)+'\t'+str(vals[0])+'\t'+str(vals[1]))


    # ## strand bias
    df_both_pass_s0 = df_both_pass_s0[ ( (df_both_pass_s0['normal_reads1_plus']>=NORMAL_STRAND_BIAS_PLUS) &
                                         (df_both_pass_s0['normal_reads1_minus']>=NORMAL_STRAND_BIAS_MINUS) &
                                         (df_both_pass_s0['tumor_reads1_plus']>=TUMOR_STRAND_BIAS_1_PLUS) &
                                         (df_both_pass_s0['tumor_reads1_minus']>=TUMOR_STRAND_BIAS_1_MINUS)  &
                                          (df_both_pass_s0['tumor_reads2_plus']>=TUMOR_STRAND_BIAS_2_PLUS) &
                                          (df_both_pass_s0['tumor_reads2_minus']>=TUMOR_STRAND_BIAS_2_MINUS)  ) ]

    print('After strand bias:'+ str(df_both_pass_s0.shape[0]))
    vals=checkRef(df_both_pass_s0,0)
    print(str(FDR)+'\t'+str(vals[0])+'\t'+str(vals[1]))


    # pval
    # sns.distplot(df_both_pass_s0['somatic_p_value'],color='green',kde=False)
    # plt.savefig(OUT_DIR+"/"+'snp_indel_somatic_pval')
    # plt.close()





    # sns.distplot(df_both_pass_s0['somatic_p_value'],color='green',kde=False)
    # plt.savefig(OUT_DIR+"/"+'snp_indel_somatic_pval_filter_3')
    # plt.close()


    # df_both_pass_s0.to_csv(OUT_DIR+"/"+SAMPLE+"_filter_final.txt",sep='\t',index=False)

    print('Filter:'+ str(df_both_pass_s0.shape[0]))
    print('Filter:snp:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='snp'].shape[0]))
    print('Filter:indel:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='indel'].shape[0]))
