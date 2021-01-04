library(ggplot2)

#change the filename below to whatever job you have
CxResults <- CIBERSORTx_Job12_Results[,0:23]

#this portion adds columns together for b cells, t cells, etc
b <- CxResults$`B cells naive` + CxResults$`B cells memory`
b_cells <- data.frame(CxResults$Mixture, b)
t_4 <- CxResults$`T cells CD4 naive` + CxResults$`T cells CD4 memory resting` +
  CxResults$`T cells CD4 memory activated`
t_cells_4 <- data.frame(CxResults$Mixture, t_4)
t_8 <- data.frame(CxResults$Mixture, CxResults$`T cells CD8`)
m <- CxResults$`Macrophages M0` + CxResults$`Macrophages M1` +
  CxResults$`Macrophages M2`
macro <- data.frame(CxResults$Mixture, m)
neutro <- data.frame(CxResults$Mixture, CxResults$Neutrophils)
d <- CxResults$`Dendritic cells activated` + CxResults$`Dendritic cells resting`
den <- data.frame(CxResults$Mixture, d)
mono <- data.frame(CxResults$Mixture, CxResults$Monocytes)
nk <- CxResults$`NK cells activated` + CxResults$`NK cells resting`
nk_cells <- data.frame(CxResults$Mixture, nk)


# bar plot for each column each mixture
# B cell, T cell CD8+, T cell CD4+, Macrophage, Neutrophil, myeloid dendritic cell
# Monocyte, NK Cell

ggplot(data=b_cells, aes(x=CxResults.Mixture, y=b)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("B Cell") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())

ggplot(data=t_cells_4, aes(x=CxResults.Mixture, y=t_4)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("T cell CD4+") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())

ggplot(data=t_8, aes(x=CxResults.Mixture, y=CxResults..T.cells.CD8.)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("T cell CD8+") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())

ggplot(data=macro, aes(x=CxResults.Mixture, y=m)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("Macrophage") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())

ggplot(data=neutro, aes(x=CxResults.Mixture, y=CxResults.Neutrophils)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("Neutrophil") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())

ggplot(data=den, aes(x=CxResults.Mixture, y=d)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("Dendritic cell") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())

ggplot(data=mono, aes(x=CxResults.Mixture, y=CxResults.Monocytes)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("Monocyte") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())

ggplot(data=nk_cells, aes(x=CxResults.Mixture, y=nk)) +
  geom_bar(stat = "identity", color = "purple", fill = "purple") +
  ggtitle("NK Cell") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.5, angle=90), axis.title.y=element_blank())




# pie plot for each sample

results <- data.frame(
  sample = CxResults$Mixture,
  B_cells = b,
  TCells_CD4 = t_4,
  TCells_CD8 = CxResults$`T cells CD8`,
  Macrophage = m,
  Neutrophil = CxResults$Neutrophils,
  Dendritic = d,
  Monocyte = CxResults$Monocytes,
  NKCell = nk,
  stringsAsFactors = FALSE
)



for (i in 1:96) {
  dafr = data.frame(
    sample = rep(results[i,1], 8),
    cell_type = c("B cells", "TCells_CD4", "TCells_CD8", "Macrophage", "Neutrophil", "Dendritic", "Monocyte", "NKCell"),
    value = c(results[i,2], results[i,3], results[i,4], results[i,5], results[i,6], results[i,7], results[i,8], results[i,9]),
    stringsAsFactors = FALSE
  )
  #print(dafr)
  plot <- ggplot(dafr, aes(x= sample, y=value, fill=cell_type)) +
            geom_bar(width = 1, stat = "identity")
  print(plot)
}

# pie plot of all samples


output <- matrix(ncol=1, nrow=768)
for(i in 1:96) {
  j = (i-1) * 8
  output[j+1, 1] <- results[i,1]
  output[j+2, 1] <- results[i,1]
  output[j+3, 1] <- results[i,1]
  output[j+4, 1] <- results[i,1]
  output[j+5, 1] <- results[i,1]
  output[j+6, 1] <- results[i,1]
  output[j+7, 1] <- results[i,1]
  output[j+8, 1] <- results[i,1]
}

samp <- data.frame(output)

c_type <-
  rep(c("B cells", "TCells_CD4", "TCells_CD8", "Macrophage", "Neutrophil", "Dendritic", "Monocyte", "NKCell"), 96)

out <- matrix(ncol = 1, nrow = 768)
for(i in 1:96) {
  y = (i-1) * 8
  out[y+1, 1] <- results[i,2]
  out[y+2, 1] <- results[i,3]
  out[y+3, 1] <- results[i,4]
  out[y+4, 1] <- results[i,5]
  out[y+5, 1] <- results[i,6]
  out[y+6, 1] <- results[i,7]
  out[y+7, 1] <- results[i,8]
  out[y+8, 1] <- results[i,9]
}

val <- data.frame(out)

df <- data.frame(
  sample = samp,
  cell_type = c_type,
  value = val,
  stringsAsFactors = FALSE
)

ggplot(df, aes(x= output, y=out, fill=cell_type)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  ggtitle("All Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 0.55, angle=90), axis.title.y=element_blank())
