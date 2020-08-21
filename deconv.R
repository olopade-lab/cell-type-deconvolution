setwd("~/Downloads/deconvGenes/") #sets directory

process <- function(file){ #function cleans data-frame by having the gene name listed as row columns rather than first column
  tpm<-read.csv(file) #reads file
  rowTpm<-as.character(tpm[,2]) 
  tpm<-tpm[,3:dim(tpm)[2]]
  row.names(tpm)<-rowTpm
  return(tpm)
}

library(immunedeconv) #loads library


allGenes<-process('./Output/allGenes.csv')

#loads processed files (tpm with hgnc symbols), then processes them
timer<-process('./Output/timer.csv')
epic<-process('./Output/epic.csv')
xcell<-process('./Output/xcell.csv')
cibersort<-process('./Output/cibersort.csv')
mcp_counter<-process('./Output/mcp_counter.csv')
cibersort_abs<-process('./Output/cibersort_abs.csv')
quantiseq<-process('./Output/quantiseq.csv')

indicVec<-rep('brca', ncol(timer)) #sets cancer type for timer algorithm

#necessary files for cibersort
set_cibersort_binary("./ciberFiles/CIBERSORT.R")
set_cibersort_mat("./ciberFiles/LM22.txt") #genes used by cibersort

#runs immune deconvolution methods 
quantiseqResults<-deconvolute(quantiseq,"quantiseq")
mcp_counterResults<-immunedeconv::deconvolute(mcp_counter,"mcp_counter")
xcellResults<-immunedeconv::deconvolute(xcell,"xcell")
epicResults<-immunedeconv::deconvolute(epic,"epic")

cibersortResults<-immunedeconv::deconvolute(cibersort,"cibersort")
cibersort_absResults<-immunedeconv::deconvolute(cibersort_abs,"cibersort_abs")

timerResults<-immunedeconv::deconvolute(timer ,"timer", indicVec)

#rows with the same information are combined (ex. different B cells are combined as 'B cell')
combine_rows <- function(data, row1, row2) {
  data[row2, 2:ncol(quantiseqResults)] <- data[row1, 2:ncol(quantiseqResults)] + data[row2, 2:ncol(quantiseqResults)]
  data[-row1, ]
}

quantiseqResults<-combine_rows(quantiseqResults,2,3)
quantiseqResults[2,1]<-"Macrophage"

cibersortResults<-combine_rows(cibersortResults,1,2)
cibersortResults<-combine_rows(cibersortResults,1,2)
cibersortResults[1,1]<-"B cell"

cibersortResults<-combine_rows(cibersortResults,3,4)
cibersortResults<-combine_rows(cibersortResults,3,4)
cibersortResults[3,1]<-"T cell CD4+"

cibersortResults<-combine_rows(cibersortResults,7,8)
cibersortResults[7,1]<-"NK cell"

cibersortResults<-combine_rows(cibersortResults,9,10)
cibersortResults<-combine_rows(cibersortResults,9,10)
cibersortResults[9,1]<-"Macrophage"

cibersortResults<-combine_rows(cibersortResults,10,11)
cibersortResults[10,1]<-"Myeloid dendritic cell"

cibersort_absResults<-combine_rows(cibersort_absResults,1,2)
cibersort_absResults<-combine_rows(cibersort_absResults,1,2)
cibersort_absResults[1,1]<-"B cell"

cibersort_absResults<-combine_rows(cibersort_absResults,3,4)
cibersort_absResults<-combine_rows(cibersort_absResults,3,4)
cibersort_absResults[3,1]<-"T cell CD4+"

cibersort_absResults<-combine_rows(cibersort_absResults,7,8)
cibersort_absResults[7,1]<-"NK cell"

cibersort_absResults<-combine_rows(cibersort_absResults,9,10)
cibersort_absResults<-combine_rows(cibersort_absResults,9,10)
cibersort_absResults[9,1]<-"Macrophage"

cibersort_absResults<-combine_rows(cibersort_absResults,10,11)
cibersort_absResults[10,1]<-"Myeloid dendritic cell"


library(ggplot2)
library(gridExtra)

#creates intersample comparison barplots

plotter<-function(df, n){
  p<-ggplot(data=df, aes(x=samples, y=values)) +
    geom_bar(stat="identity", color="blue", fill="blue") +
    theme_classic() +
    #ggtitle(n) +
    theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())
  p
}



quantiseqBcell<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[1,2:ncol(quantiseqResults)]))
quantiseqTcell8<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[7,2:ncol(quantiseqResults)]))
quantiseqTcell4<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[6,2:ncol(quantiseqResults)]))
quantiseqMacrophage<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[2,2:ncol(quantiseqResults)]))
quantiseqNeutrophil<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[4,2:ncol(quantiseqResults)]))
quantiseqmDC<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[9,2:ncol(quantiseqResults)]))
quantiseqMonocyte<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[3,2:ncol(quantiseqResults)]))
quantiseqNKcell<-data.frame(samples=c(colnames(quantiseqResults)[2:ncol(quantiseqResults)]), values=as.numeric(quantiseqResults[5,2:ncol(quantiseqResults)]))

quantiseqBcell<-plotter(quantiseqBcell,"B Cell")
quantiseqTcell8<-plotter(quantiseqTcell8, "T cell CD8+")
quantiseqTcell4<-plotter(quantiseqTcell4, "T cell CD4+")
quantiseqMacrophage<-plotter(quantiseqMacrophage, "Macrophage")
quantiseqNeutrophil<-plotter(quantiseqNeutrophil, "Neutrophil")
quantiseqmDC<-plotter(quantiseqmDC, "Myeloid dendritic cell")
quantiseqMonocyte<-plotter(quantiseqMonocyte, "Monocyte")
quantiseqNKcell<-plotter(quantiseqNKcell, "NK cell")



timerBcell<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=as.numeric(timerResults[1,2:ncol(timerResults)]))
timerTcell8<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=as.numeric(timerResults[3,2:ncol(timerResults)]))
timerTcell4<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=as.numeric(timerResults[2,2:ncol(timerResults)]))
timerMacrophage<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=as.numeric(timerResults[5,2:ncol(timerResults)]))
timerNeutrophil<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=as.numeric(timerResults[4,2:ncol(timerResults)]))
timermDC<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=as.numeric(timerResults[6,2:ncol(timerResults)]))
timerMonocyte<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=rep(0, ncol(timerResults)-1))
timerNKcell<-data.frame(samples=c(colnames(timerResults)[2:ncol(timerResults)]), values=rep(0, ncol(timerResults)-1))

timerBcell<-plotter(timerBcell,"B Cell")
timerTcell8<-plotter(timerTcell8, "T cell CD8+")
timerTcell4<-plotter(timerTcell4, "T cell CD4+")
timerMacrophage<-plotter(timerMacrophage, "Macrophage")
timerNeutrophil<-plotter(timerNeutrophil, "Neutrophil")
timermDC<-plotter(timermDC, "Myeloid dendritic cell")
timerMonocyte<-plotter(timerMonocyte, "Monocyte")
timerNKcell<-plotter(timerNKcell, "NK cell")

#cibersort

cibersortBcell<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[1,2:ncol(cibersortResults)]))
cibersortTcell8<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[2,2:ncol(cibersortResults)]))
cibersortTcell4<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[3,2:ncol(cibersortResults)]))
cibersortMacrophage<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[9,2:ncol(cibersortResults)]))
cibersortNeutrophil<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[14,2:ncol(cibersortResults)]))
cibersortmDC<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[10,2:ncol(cibersortResults)]))
cibersortMonocyte<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[8,2:ncol(cibersortResults)]))
cibersortNKcell<-data.frame(samples=c(colnames(cibersortResults)[2:ncol(cibersortResults)]), values=as.numeric(cibersortResults[7,2:ncol(cibersortResults)]))

cibersortBcell<-plotter(cibersortBcell,"B Cell")
cibersortTcell8<-plotter(cibersortTcell8, "T cell CD8+")
cibersortTcell4<-plotter(cibersortTcell4, "T cell CD4+")
cibersortMacrophage<-plotter(cibersortMacrophage, "Macrophage")
cibersortNeutrophil<-plotter(cibersortNeutrophil, "Neutrophil")
cibersortmDC<-plotter(cibersortmDC, "Myeloid dendritic cell")
cibersortMonocyte<-plotter(cibersortMonocyte, "Monocyte")
cibersortNKcell<-plotter(cibersortNKcell, "NK cell")

#cibersort_abs

cibersort_absBcell<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[1,2:ncol(cibersort_absResults)]))
cibersort_absTcell8<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[2,2:ncol(cibersort_absResults)]))
cibersort_absTcell4<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[3,2:ncol(cibersort_absResults)]))
cibersort_absMacrophage<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[9,2:ncol(cibersort_absResults)]))
cibersort_absNeutrophil<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[14,2:ncol(cibersort_absResults)]))
cibersort_absmDC<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[10,2:ncol(cibersort_absResults)]))
cibersort_absMonocyte<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[8,2:ncol(cibersort_absResults)]))
cibersort_absNKcell<-data.frame(samples=c(colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]), values=as.numeric(cibersort_absResults[7,2:ncol(cibersort_absResults)]))

cibersort_absBcell<-plotter(cibersort_absBcell,"B Cell")
cibersort_absTcell8<-plotter(cibersort_absTcell8, "T cell CD8+")
cibersort_absTcell4<-plotter(cibersort_absTcell4, "T cell CD4+")
cibersort_absMacrophage<-plotter(cibersort_absMacrophage, "Macrophage")
cibersort_absNeutrophil<-plotter(cibersort_absNeutrophil, "Neutrophil")
cibersort_absmDC<-plotter(cibersort_absmDC, "Myeloid dendritic cell")
cibersort_absMonocyte<-plotter(cibersort_absMonocyte, "Monocyte")
cibersort_absNKcell<-plotter(cibersort_absNKcell, "NK cell")

#mcp_counter

mcp_counterBcell<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=as.numeric(mcp_counterResults[5,2:ncol(mcp_counterResults)]))
mcp_counterTcell8<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=as.numeric(mcp_counterResults[2,2:ncol(mcp_counterResults)]))
mcp_counterTcell4<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=rep(0, ncol(mcp_counterResults)-1))
mcp_counterMacrophage<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=as.numeric(mcp_counterResults[7,2:ncol(mcp_counterResults)]))
mcp_counterNeutrophil<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=as.numeric(mcp_counterResults[9,2:ncol(mcp_counterResults)]))
mcp_countermDC<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=as.numeric(mcp_counterResults[8,2:ncol(mcp_counterResults)]))
mcp_counterMonocyte<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=as.numeric(mcp_counterResults[6,2:ncol(mcp_counterResults)]))
mcp_counterNKcell<-data.frame(samples=c(colnames(mcp_counterResults)[2:ncol(mcp_counterResults)]), values=as.numeric(mcp_counterResults[4,2:ncol(mcp_counterResults)]))

mcp_counterBcell<-plotter(mcp_counterBcell,"B Cell")
mcp_counterTcell8<-plotter(mcp_counterTcell8, "T cell CD8+")
mcp_counterTcell4<-plotter(mcp_counterTcell4, "T cell CD4+")
mcp_counterMacrophage<-plotter(mcp_counterMacrophage, "Macrophage")
mcp_counterNeutrophil<-plotter(mcp_counterNeutrophil, "Neutrophil")
mcp_countermDC<-plotter(mcp_countermDC, "Myeloid dendritic cell")
mcp_counterMonocyte<-plotter(mcp_counterMonocyte, "Monocyte")
mcp_counterNKcell<-plotter(mcp_counterNKcell, "NK cell")

#xcell


xcellBcell<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[2,2:ncol(xcellResults)]))
xcellTcell8<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[9,2:ncol(xcellResults)]))
xcellTcell4<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[3,2:ncol(xcellResults)]))
xcellMacrophage<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[21,2:ncol(xcellResults)]))
xcellNeutrophil<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[28,2:ncol(xcellResults)]))
xcellmDC<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[15,2:ncol(xcellResults)]))
xcellMonocyte<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[26,2:ncol(xcellResults)]))
xcellNKcell<-data.frame(samples=c(colnames(xcellResults)[2:ncol(xcellResults)]), values=as.numeric(xcellResults[29,2:ncol(xcellResults)]))

xcellBcell<-plotter(xcellBcell,"B Cell")
xcellTcell8<-plotter(xcellTcell8, "T cell CD8+")
xcellTcell4<-plotter(xcellTcell4, "T cell CD4+")
xcellMacrophage<-plotter(xcellMacrophage, "Macrophage")
xcellNeutrophil<-plotter(xcellNeutrophil, "Neutrophil")
xcellmDC<-plotter(xcellmDC, "Myeloid dendritic cell")
xcellMonocyte<-plotter(xcellMonocyte, "Monocyte")
xcellNKcell<-plotter(xcellNKcell, "NK cell")

#epic

epicBcell<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=as.numeric(epicResults[1,2:ncol(epicResults)]))
epicTcell8<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=as.numeric(epicResults[4,2:ncol(epicResults)]))
epicTcell4<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=as.numeric(epicResults[3,2:ncol(epicResults)]))
epicMacrophage<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=as.numeric(epicResults[6,2:ncol(epicResults)]))
epicNeutrophil<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=rep(0, ncol(epicResults)-1))
epicmDC<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=rep(0, ncol(epicResults)-1))
epicMonocyte<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=rep(0, ncol(epicResults)-1))
epicNKcell<-data.frame(samples=c(colnames(epicResults)[2:ncol(epicResults)]), values=as.numeric(epicResults[7,2:ncol(epicResults)]))

epicBcell<-plotter(epicBcell,"B Cell")
epicTcell8<-plotter(epicTcell8, "T cell CD8+")
epicTcell4<-plotter(epicTcell4, "T cell CD4+")
epicMacrophage<-plotter(epicMacrophage, "Macrophage")
epicNeutrophil<-plotter(epicNeutrophil, "Neutrophil")
epicmDC<-plotter(epicmDC, "Myeloid dendritic cell")
epicMonocyte<-plotter(epicMonocyte, "Monocyte")
epicNKcell<-plotter(epicNKcell, "NK cell")

g<-arrangeGrob(arrangeGrob(arrangeGrob(quantiseqBcell, top="B cell"), arrangeGrob(quantiseqTcell8, top="T Cell CD8+"), arrangeGrob(quantiseqTcell4, top="T Cell CD4+"), arrangeGrob(quantiseqMacrophage, top="Macrophage"), arrangeGrob(quantiseqNeutrophil, top="Neutrophil"), arrangeGrob(quantiseqmDC, top="Myeloid dendritic cell"), arrangeGrob(quantiseqMonocyte, top="Monocyte"), arrangeGrob(quantiseqNKcell, top="NK Cell"), left="quanTIseq", nrow=1),
               arrangeGrob(timerBcell, timerTcell8, timerTcell4, timerMacrophage, timerNeutrophil, timermDC, timerMonocyte, timerNKcell, left="timer", nrow=1), 
               arrangeGrob(cibersortBcell, cibersortTcell8, cibersortTcell4, cibersortMacrophage, cibersortNeutrophil, cibersortmDC, cibersortMonocyte, cibersortNKcell, left="cibersort", nrow=1),
               arrangeGrob(cibersort_absBcell, cibersort_absTcell8, cibersort_absTcell4, cibersort_absMacrophage, cibersort_absNeutrophil, cibersort_absmDC, cibersort_absMonocyte, cibersort_absNKcell, left="cibersort_abs", nrow=1),
               arrangeGrob(mcp_counterBcell, mcp_counterTcell8, mcp_counterTcell4, mcp_counterMacrophage, mcp_counterNeutrophil, mcp_countermDC, mcp_counterMonocyte, mcp_counterNKcell, left="mcp_counter", nrow=1),
               arrangeGrob(xcellBcell, xcellTcell8, xcellTcell4, xcellMacrophage, xcellNeutrophil, xcellmDC, xcellMonocyte, xcellNKcell, left="xcell", nrow=1),
               arrangeGrob(epicBcell, epicTcell8, epicTcell4, epicMacrophage, epicNeutrophil, epicmDC, epicMonocyte, epicNKcell, left="epic", nrow=1),
               nrow=7)
ggsave('./plots/barplot.pdf',g, height=15, width = 50, limitsize = FALSE)


library(dplyr)
library(tidyr)
plottingFunctions<-function(results, name){
  
  #for (i in colnames(results)[2:84]) {
  #name<-paste(i, deparse(substitute(results)), sep = "_")
  #assign(name, results[,c("cell_type",i)])
  #assign("res", results[,c("cell_type",i)])
  
  results %>%
    gather(sample, fraction, -cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    ylab('quantity') +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(results)))
  #  fileN=paste("./deconvFiguresSamples/",name,sep = "_")
  #  fileN=paste(fileN,"png",sep = ".")
  #  ggsave(fileN)
  ggsave(paste('./plots/',name,'.pdf', sep=""), height=15, width = 15, limitsize = FALSE)
  #}
}

plottingFunctions(epicResults, 'epic')
plottingFunctions(quantiseqResults, 'quantiseq')
plottingFunctions(mcp_counterResults,'mcp_counter')
#plottingFunctions(xcellResults,'xcell' )
#plottingFunctions(ciberResults, 'ciber')
#plottingFunctions(ciber_absResults, 'ciber_abs')
plottingFunctions(timerResults, 'timer')





library(ggplot2)
library(scales)
library(gridExtra)
#creates intrasample comparison bar plots

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text = element_blank(),
    panel.grid  = element_blank(),
    legend.title = element_blank(),
    #legend.position = "none"
  )

pieplot<-function(df, name){
  
  
  
  bp<- ggplot(df, aes(x="", y=value, fill=cell_type)) +
    geom_bar(width = 1, stat = "identity") +
    ggtitle(name) +
    blank_theme
  bp
  #  pie <- bp + coord_polar("y", start=0) +
  #    blank_theme
  
  #  pie
}

myplotsCibersort <- vector('list', 0)
myplotsCibersort_abs <- vector('list', 0)
myplotsQuantiseq <- vector('list', 0)
myplotsepic <- vector('list', 0)

for (i in colnames(cibersortResults)[2:ncol(cibersortResults)]) {
  df<-data.frame(cell_type=cibersortResults$cell_type, value=cibersortResults[[i]])
  
  myplotsCibersort[[i]] <- local({
    i <- i
    p<-pieplot(df, i)
    print(p)
  })
  
}

for (i in colnames(cibersort_absResults)[2:ncol(cibersort_absResults)]) {
  df<-data.frame(cell_type=cibersort_absResults$cell_type, value=cibersort_absResults[[i]])
  
  myplotsCibersort_abs[[i]] <- local({
    i <- i
    p<-pieplot(df, i)
    print(p)
  })
  
}

for (i in colnames(quantiseqResults)[2:ncol(quantiseqResults)]) {
  df<-data.frame(cell_type=quantiseqResults$cell_type, value=quantiseqResults[[i]])
  
  myplotsQuantiseq[[i]] <- local({
    i <- i
    p<-pieplot(df, i)
    print(p)
  })
  
}

for (i in colnames(epicResults)[2:ncol(epicResults)]) {
  df<-data.frame(cell_type=epicResults$cell_type, value=epicResults[[i]])
  
  myplotsepic[[i]] <- local({
    i <- i
    p<-pieplot(df, i)
    print(p)
  })
  
}


g<-arrangeGrob(grobs = c(myplotsCibersort, myplotsCibersort_abs, myplotsQuantiseq, myplotsepic), nrow = 4)
ggsave('./plots/pieplot.pdf',g, height=20, width = 500, limitsize = FALSE)



