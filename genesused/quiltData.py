#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 17:16:19 2020

@author: arvindkrishnan
"""

import pandas as pd

#list1=pd.read_csv('/Users/arvindkrishnan/Downloads/deconvGenes/xcell.txt', header=None) #xcell
list1=pd.read_csv('./xcellR.txt')
list1=list1['x']
list1=list1.to_frame()


list2=pd.read_csv('./mcpCounter.txt', sep="\t") #mcpCounter
list2=list2['HUGO symbols']
list2=list2.to_frame()

list3=pd.read_csv('./Cibersort.txt', sep="\t") #cibersort
list3=list3['Gene symbol']
list3=list3.to_frame()

list4=pd.read_csv('./timer.txt', sep="\t") #timer

list5=pd.read_csv('./quantiseq.txt', sep="\t") #quantiseq
list5=list5['ID']
list5=list5.to_frame()

list6=pd.read_csv('./epic_BRef.txt', sep="\t")#epic
list7=pd.read_csv('./epic_TRef.txt', sep="\t")






list1=list1.rename(columns={'x':'gene'})
listA=list1
#listA=listA.rename(columns={0:'gene'})
list2=list2.rename(columns={'HUGO symbols':'gene'})
list3=list3.rename(columns={'Gene symbol':'gene'})
list4=list4.rename(columns={'x':'gene'})
list5=list5.rename(columns={'ID':'gene'})
list6=list6.rename(columns={'x':'gene'})
list7=list7.rename(columns={'x':'gene'})

list1.to_csv('./xcell.csv')
list2.to_csv('./mcpCounter.csv')
list3.to_csv('./Cibersort.csv')
list4.to_csv('./timer.csv')
list5.to_csv('./quantiseq.csv')

listEpic=list6
listEpic=listEpic.append(list7)
listEpic.to_csv('./epic.csv')

listA=list1
listA=listA.append(list2)
listA=listA.append(list3)
listA=listA.append(list4)
listA=listA.append(list5)
listA=listA.append(list6)
listA=listA.append(list7)
listA.to_csv('./allGenesUsed.csv')
#listA=listA['gene'].to_list()



