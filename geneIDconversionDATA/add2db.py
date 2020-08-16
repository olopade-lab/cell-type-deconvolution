#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 18:20:17 2020

@author: arvindkrishnan
"""
import pandas as pd

database=pd.read_csv('./database.csv')

def AppendDatabase(df, database):
    df=df[['InputTerm','Ensembl']]
    df2=pd.DataFrame({'Approved symbol':df['InputTerm'],'Ensembl gene ID': df['Ensembl']})
    database=database.append(df2)
    return(database)

df=pd.read_excel('./GeneALaCart-4142-200720-214249.xlsx%3Fsv=2015-04-05&sr=b&sig=noEmyoxZCJHy3qlBauPW025iarMYHuSV5MAT57EKUuE=&st=2020-07-20T22%3A51%3A43Z&se=2020-07-21T02%3A56%3A43Z&sp=r.xlsx', sheet_name='ExternalIdentifiers')
database=AppendDatabase(df, database)

df=pd.read_excel('./GeneALaCart-500-200811-010737.xlsx%3Fsv=2015-04-05&sr=b&sig=7FZrlH7Dn9H3PHA4KUF0Q+gA6R4gH1Va8n9d6xu9tuU=&st=2020-08-11T01%3A07%3A19Z&se=2020-08-11T05%3A12%3A19Z&sp=r.xlsx', sheet_name='ExternalIdentifiers')
database=AppendDatabase(df, database)

df=pd.read_excel('./GeneALaCart-500-200811-011258.xlsx%3Fsv=2015-04-05&sr=b&sig=H+BDmC5WfAAYgyIMDwiDnQUTo7DEqPp9f0NiNZpxQeU=&st=2020-08-11T01%3A10%3A01Z&se=2020-08-11T05%3A15%3A01Z&sp=r.xlsx', sheet_name='ExternalIdentifiers')
database=AppendDatabase(df, database)

df=pd.read_excel('./GeneALaCart-500-200811-011738.xlsx%3Fsv=2015-04-05&sr=b&sig=sRCEFBa%2Ffkah5vciGvMWFAykGhDO8FZwNNAUEQVGf3o=&st=2020-08-11T01%3A13%3A57Z&se=2020-08-11T05%3A18%3A57Z&sp=r.xlsx', sheet_name='ExternalIdentifiers')
database=AppendDatabase(df, database)

database.to_csv('./database.csv')