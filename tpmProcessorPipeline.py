#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 12:59:20 2020

@author: arvindkrishnan
"""

import argparse 


import pandas as pd
import numpy as np

def converter(ensg, database, listA):
    '''
    

    Parameters
    ----------
    ensg : df
        column of ensembl ids from tpm.
    database : df
        df mapping ensembl to hgnc ids.
    listA : list
        list of genes used.

    Returns
    -------
    newVal : 
        hgnc ids.
    place : 
        Whether this is an alias or not

    '''
    sym=[]
    place='none'
    sym1=database[database['Ensembl gene ID']==ensg]['Approved symbol'] #finds hgnc symbol from ensembl ID
    if len(sym1)>0: #if such hgnc id exists
        sym1=sym1[sym1.index[0]] #extracts the symbol from the row
        sym.append(sym1) #stores this symbol
    else:
        sym1='none'
    sym2=database[database['Ensembl gene ID']==ensg]['Alias symbol'] #finds alias hgnc symbol
    if len(sym2)>0:
        for i in range(0, len(sym2)):
            sym.append(sym2[sym2.index[i]]) #appends all aliases
    else:
        sym2='none'
        
    if len(sym)>0:
        for c, i in enumerate(sym):
            if i in listA: #this loop includes the hgnc id that is recognized by the algorithm
                newVal=i
                if c==0:
                    place='approved'
                else:
                    place='alias'
                break
            else:
                newVal=np.nan
    else:
        newVal=np.nan
    return newVal, place


def BACKconverter(hgnc, database):
    '''
    
    Function isn't used but operates similar to converter (above). 
    Changes hgnc to ensembl

    '''
    sym=[]
    sym1=database[database['Approved symbol']==hgnc]['Ensembl gene ID']
    if len(sym1)>0:
        sym1=sym1[sym1.index[0]]
        sym.append(sym1)
    else:
        sym1='none'
    sym2=database[database['Alias symbol']==hgnc]['Ensembl gene ID']
    if len(sym2)>0:
        for i in range(0, len(sym2)):
            sym.append(sym2[sym2.index[i]])
    else:
        sym2='none'

    if len(sym)>0:
        newVal=sym[0]
    else:
        newVal=np.nan
    return newVal


#tpmensg
def dup(tpm, tpmensg):
    duplicates=tpm[tpm['gene_id'].duplicated(keep=False)==True]
    duplicates=duplicates.sort_values(by=['gene_id'])
    duplicates=duplicates[['gene_id','where']]
    #duplicates=duplicates.to_frame()
    duplicates['ensg']=tpmensg[tpmensg.index.isin(duplicates.index)]
    return duplicates






def main():
    '''
    main function
    
    Command Line arguements: 
        tpm=tpm nomalized expression file
        method=algorithm to use
        -o outputname
        
    

    Returns
    -------
    None.

    '''
    #Command line arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('tpm', help='tpm file')
    parser.add_argument('method', help='algorithm to be used: all, quantiseq, mcp_counter, epic...')
    parser.add_argument('-o','--outputname', default= 'TpmProcessed', help='name of output file')
    args = parser.parse_args()
    
    #Loads csv of genes used by algorithms depending on method selected
    if args.method=='all':
        listA=pd.read_csv('./genesused/allGenesUsed.csv')
    elif args.method=='quantiseq':
        listA=pd.read_csv('./genesused/quantiseq.csv')
    elif args.method=='mcp_counter':
        listA=pd.read_csv('./genesused/mcpCounter.csv')
    elif args.method=='cibersort':
        listA=pd.read_csv('./genesused/Cibersort.csv')
    elif args.method=='cibersort_abs':
        listA=pd.read_csv('./genesused/Cibersort.csv')
    elif args.method=='timer':
        listA=pd.read_csv('./genesused/timer.csv')
    elif args.method=='xcell':
        listA=pd.read_csv('./genesused/xcell.csv')
    elif args.method=='epic':
        listA=pd.read_csv('./genesused/epic.csv')
        
    #changes dataframe to list
    listA=listA['gene'].to_list()

    #Loads dataframe that has ennsembl ids, corresponding hngc gene symbols, and aliases. 
    database=pd.read_csv('./geneIDconversionDATA/database.csv')

    #reads tpm normalized expression profile inputted as command line arguement
    tpm=pd.read_csv(args.tpm)
    #Changes ensemblIDs to hgnc gene symbols
    a=tpm.apply(lambda x: converter(x['gene_id'], database, listA ), axis=1)
    tpm['gene_id'], tpm['where']=zip(*a)

    #reads genes in tpm file
    tpmensg=pd.read_csv(args.tpm)['gene_id']
    
    #drops na values, i.e genes not used by algorithms
    tpm=tpm.dropna()
    
    #Finds duplicates
    duplicates=dup(tpm, tpmensg)

    #to remove
    tpm=tpm.drop(duplicates[duplicates['where']=='alias'].index)
    duplicates=dup(tpm, tpmensg)
    tpm=tpm.drop(['where'],axis=1)

    #duplicates are written out as well as output file
    duplicates.to_csv('./Output/duplicates.csv')
    tpm.to_csv('./Output/'+ args.outputname + '.csv', index_label='gene_id')
    
if __name__ == "__main__":
    main()
    