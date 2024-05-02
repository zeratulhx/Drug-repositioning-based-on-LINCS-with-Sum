import pandas as pd
from tqdm import tqdm
from cmapPy.pandasGEXpress.parse import parse
import argparse
from rapidfuzz import process, fuzz
parser1=argparse.ArgumentParser(description='Specify file paths')
parser1.add_argument('--sig_info',type=str,default='siginfo_beta.txt')
parser1.add_argument('--gene_info',type=str,default='geneinfo_beta.txt')
parser1.add_argument('--query',type=str,default='query.xlsx')
parser1.add_argument('--batch_size',type=int,default=1000)
parser1.add_argument('--trt_cp',type=str,default='level5_beta_trt_cp_n720216x12328.gctx')
parser1.add_argument('--Xsum_num',type=int,default=200)
parser1.add_argument('--cmp_info',type=str,default='compoundinfo_beta.txt')
parser1.add_argument('--drug_info',type=str,default='repurposing_drugs_20200324.txt')
parser1.add_argument('--output',type=str,default='repositioned_drug.csv')
parser1.add_argument('--cell',type=str,default='AG06263_2')



args=parser1.parse_args()
#This function process the filter between perturbations.
def get_cid(sig_info):
    #compounds only
    sig_info_df=sig_info[sig_info['pert_type']=='trt_cp']
    #concentration
    sig_info_df=sig_info_df[sig_info_df['pert_dose_unit']=='uM']
    sig_info_df=sig_info_df[sig_info_df['nearest_dose']==10]
    #perturbation time
    sig_info_df=sig_info_df[sig_info_df['pert_time']==6]
    #cell line
    sig_info_df=sig_info_df[sig_info_df['cell_iname']==args.cell]
    #high quality signatures only
    sig_info_df=sig_info_df[sig_info_df['is_hiq']==1]
    return(sig_info_df['sig_id'])


#This function extracts the landmark genes and overlapping genes from query list for gctx parsing.
def get_rid(geneinfo_df,query):
    gene_overlapping=geneinfo_df[geneinfo_df['gene_symbol'].isin(query['symbol'])]
    gene_overlapping=gene_overlapping[gene_overlapping['feature_space']=='landmark']
    return(gene_overlapping['gene_id'].map(str))





#Thisd function extracts top N expression profiles and calculate XSum score.
def Xsum(col,N,geneinfo_df,query):
    pos=col.head(N)
    neg=col.tail(N)
    changed=pd.concat([pos,neg],axis=0)
    num=[0,0]
    for i in range(0,2):
        term=['Up','Down'][i]
        valid_indices=[]
        gene_info_paper=geneinfo_df[geneinfo_df['gene_symbol'].isin(query[query['orientation']==term]['symbol'])]
        #print(gene_info_paper.head(200)['gene_id'])
        for index in gene_info_paper.head(200)['gene_id']:
            index1=str(index)
            if(index1 in changed.index):
                valid_indices.append(index1)
        num[i]=sum(changed[valid_indices])

    return(num[0]-num[1])
    
def main():
    try:
        sig_info=pd.read_csv(args.sig_info,sep='\t')
    except:
        print('Invalid siginfo file.')
    cid=get_cid(sig_info)
    try:
        print(args.gene_info)
        geneinfo_df=pd.read_csv(args.gene_info,sep='\t')
    except:
        print('Invalid geneinfo file.')
    
    try:
        query=pd.read_excel(args.query)
    except:
        print('Invalid query file.')
    
    rids=get_rid(geneinfo_df,query)
    #This block parses the gctx file and generate sorted XSum results
    df_output=pd.DataFrame()
    batch_size=args.batch_size
    for start_idx in tqdm(range(0,len(cid),batch_size)):
        end_idx=min(start_idx+batch_size,len(cid))
        batch_sig_ids=cid[start_idx:end_idx]

        pert_data=parse(args.trt_cp,cid=batch_sig_ids,rid=rids)
        for col in pert_data.data_df.columns:
            sorted_col=pert_data.data_df[col].sort_values(ascending=False)
            score=Xsum(sorted_col,args.Xsum_num,geneinfo_df,query)
            line=len(df_output)
            df_output.loc[line,'pert']=col
            df_output.loc[line,'score']=score
    df_output.sort_values(by='score',ascending=True,inplace=True)
    del(pert_data)

    #This block extracts compound info data(names moa) for further annotation
    df_comp=pd.read_csv(args.cmp_info,sep='\t')
    mapping = sig_info.set_index('sig_id')['pert_id'].to_dict()
    df_output['pert_name'] = df_output['pert'].map(mapping)
    mapping = df_comp.set_index('pert_id')['compound_aliases'].to_dict()
    df_output['aliases'] = df_output['pert_name'].map(mapping)
    mapping = df_comp.set_index('pert_id')['cmap_name'].to_dict()
    df_output['name'] = df_output['pert_name'].map(mapping)
    mapping = df_comp.set_index('pert_id')['moa'].to_dict()
    df_output['moa'] = df_output['pert_name'].map(mapping)
    df_drug=pd.read_csv(args.drug_info,sep='\t',skiprows=9)
    df_drug['pert_iname']=df_drug['pert_iname'].str.lower()

    #This block annotates compounds with clinical phase data
    for index,row in df_output.iterrows():
        try:
            name1=row['name']
        except:
            name1=''
        try:
            name2=row['aliases']
        except:
            name2=''
        if name1 in df_drug['pert_iname'].values.tolist() :
            df_output.loc[index,'phase']=df_drug[df_drug['pert_iname']==name1].iloc[0]['clinical_phase']
        elif name2 in df_drug['pert_iname'].values.tolist():
            df_output.loc[index,'phase']=df_drug[df_drug['pert_iname']==name2].iloc[0]['clinical_phase']
        else:
            try:
                pair1, score1,_ = process.extractOne(name1, df_drug['pert_iname'].values.tolist(), scorer=fuzz.WRatio)
            except:
                score1=0
            try:
                pair2, score2,_ = process.extractOne(name2, df_drug['pert_iname'].values.tolist(), scorer=fuzz.WRatio)
            except:
                score2=0
            if score1 > 90 or score2 > 90:  
                if score1 > score2:
                    matched_name = pair1
                    score = score1
                else:
                    matched_name = pair2
                    score = score2
                if(score1 ==100 or score2==100 ):
                    df_output.loc[index, 'is_poss_pair'] = 1

                df_output.loc[index, 'phase'] = df_drug[df_drug['pert_iname'] == matched_name].iloc[0]['clinical_phase']
                df_output.loc[index, 'pair'] = matched_name
                df_output.loc[index, 'pair_score'] = score
    df_output.to_csv(args.output,sep=',',index=False)


    

if __name__=='__main__':
    main()
