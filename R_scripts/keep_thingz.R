
### PCOA


ordi23s <- ordinate(file2meco::meco2phyloseq(data_23s_hell), "PCoA", "horn")

expl1 <- ordi23s$values$Relative_eig[1]
expl2 <- ordi23s$values$Relative_eig[2]
expl3 <- ordi23s$values$Relative_eig[3]
df_ordi_23s <- data.frame(dim1=ordi23s$vectors[,1],
                          dim2=ordi23s$vectors[,2],
                          dim3=ordi23s$vectors[,3],
                          country=data_23s_hell$sample_table$Country
)


pcoa23s <- ggplot(df_ordi_23s,aes(x=dim1,y=dim2,color=country))+
    geom_point()+
    theme_classic()+
    ylab(paste0("PCoA 2 (",round(expl2,digits = 3)*100,"%)"))+
    xlab(paste0("PCoA 1 (",round(expl1,digits = 3)*100,"%)"))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    ggtitle("MAPP_23s")+
    theme(plot.title = element_text(face="bold"))

ordi16s <- ordinate(file2meco::meco2phyloseq(data_16s_hell), "PCoA", "horn")

expl1 <- ordi16s$values$Relative_eig[1]
expl2 <- ordi16s$values$Relative_eig[2]
expl3 <- ordi16s$values$Relative_eig[3]
df_ordi_16s <- data.frame(dim1=ordi16s$vectors[,1],
                          dim2=ordi16s$vectors[,2],
                          dim3=ordi16s$vectors[,3],
                          country=data_16s_hell$sample_table$Country
)


pcoa16s <- ggplot(df_ordi_16s,aes(x=dim1,y=dim2,color=country))+
    geom_point()+
    theme_classic()+
    ylab(paste0("PCoA 2 (",round(expl2,digits = 3)*100,"%)"))+
    xlab(paste0("PCoA 1 (",round(expl1,digits = 3)*100,"%)"))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    ggtitle("MAPP_16s")+
    theme(plot.title = element_text(face="bold"))


ordi18s <- ordinate(file2meco::meco2phyloseq(data_18s_hell), "PCoA", "horn")

expl1 <- ordi18s$values$Relative_eig[1]
expl2 <- ordi18s$values$Relative_eig[2]
expl3 <- ordi18s$values$Relative_eig[3]
df_ordi_18s <- data.frame(dim1=ordi18s$vectors[,1],
                          dim2=ordi18s$vectors[,2],
                          dim3=ordi18s$vectors[,3],
                          country=data_18s_hell$sample_table$Country
)


pcoa18s <- ggplot(df_ordi_18s,aes(x=dim1,y=dim2,color=country))+
    geom_point()+
    theme_classic()+
    ylab(paste0("PCoA 2 (",round(expl2,digits = 3)*100,"%)"))+
    xlab(paste0("PCoA 1 (",round(expl1,digits = 3)*100,"%)"))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    ggtitle("MAPP_18s")+
    theme(plot.title = element_text(face="bold"))
