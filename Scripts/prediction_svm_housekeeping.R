#R script to predict gene expression based on Support Vector Machines,
#
# SVM are computed by e1071 package
#
# version: 0.1 06.12.2016: Enhanced plots, adapt plot limits to data (99th percentile)


library(e1071)
args <- commandArgs()
data_broad<-read.csv(args[4],header=TRUE,sep="\t",stringsAsFactors=FALSE)
data_small<-read.csv(args[5],header=TRUE,sep="\t",stringsAsFactors=FALSE)
output_dir<-args[6]
ref_dir<-args[7]
data<-data.frame(Gene=data_broad$Gene,TSS=data_broad$TSS,Broad=data_broad$Coverage,Small=data_small$Coverage)
#############################################################################################################################################################################
analyze_and_plot_clustering<-function(data,output_dir) {
    cols<-rep("black",length(data$Type))
    cols[which(data$Clustering == "Expressed")]<-"red"
    cols[which(data$Clustering == "Unexpressed")]<-"green"
    png(paste(output_dir,"/Broad_vs_Small_TSS_coverage_clustering.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,30),ylim=c(0,30))
    dev.off()

    for (xx in which(data$Clustering == "Unexpressed")) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in which(data$Clustering == "Expressed")) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

}
#############################################################################################################################################################################
sensitivity_top5000<-function(data,output_dir){
    #Check in more detail whether Top5000 genes are in expressde cluster and bottom 1000 in unexpressed cluster:
    top5000<-read.csv(paste(ref_dir,"/Top5000_NMonly.txt",sep=""),header=FALSE)
    bottom5000<-read.csv(paste(ref_dir,"/Bottom5000_NMonly.txt",sep=""),header=FALSE)
    broad_boundary<-quantile(data$Broad,0.99)
    small_boundary<-quantile(data$Small,0.99)
    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top5000$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top5000"
    }
    for (bottom in bottom5000$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom5000"
    }
    true_expressed<-which(expressed == "Top5000" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom5000" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom5000" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top5000" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0.8,0.8,0.8,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(215,25,28,maxColorValue=255)
    cols[true_unexpressed]<-rgb(44,123,182,maxColorValue=255)
    cols[false_expressed]<-rgb(171,217,233,maxColorValue=255)
    cols[false_unexpressed]<-rgb(253,174,97,maxColorValue=255)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top5000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top5000_.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()


    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top5000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top5000_Highlight.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    points(data$Broad[true_expressed],data$Small[true_expressed],col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],col=cols[false_expressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top5000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom5000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top5000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom5000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    write.table(t(c("Top5000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom5000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top5000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom5000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}

#############################################################################################################################################################################
sensitivity_top1000<-function(data,output_dir){
    #Check in more detail whether Top1000 genes are in expressde cluster and bottom 1000 in unexpressed cluster:
    top1000<-read.csv(paste(ref_dir,"/Top1000_NMonly.txt",sep=""),header=FALSE)
    bottom1000<-read.csv(paste(ref_dir,"/Bottom1000_NMonly.txt",sep=""),header=FALSE)


    broad_boundary<-quantile(data$Broad,0.99)
    small_boundary<-quantile(data$Small,0.99)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top1000$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top1000"
    }
    for (bottom in bottom1000$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom1000"
    }
    true_expressed<-which(expressed == "Top1000" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom1000" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom1000" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top1000" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0.8,0.8,0.8,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(215,25,28,maxColorValue=255)
    cols[true_unexpressed]<-rgb(44,123,182,maxColorValue=255)
    cols[false_expressed]<-rgb(171,217,233,maxColorValue=255)
    cols[false_unexpressed]<-rgb(253,174,97,maxColorValue=255)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(expressed == "Top1000" & data$Clustering == "Unexpressed")]<-"red"
    cols[which(expressed == "Bottom1000" & data$Clustering == "Expressed")]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()

    cols<-rep(rgb(0.8,0.8,0.8,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(215,25,28,maxColorValue=255)
    cols[true_unexpressed]<-rgb(44,123,182,maxColorValue=255)
    cols[false_expressed]<-rgb(171,217,233,maxColorValue=255)
    cols[false_unexpressed]<-rgb(253,174,97,maxColorValue=255)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()

    cols<-rep(rgb(0.8,0.8,0.8,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(215,25,28,maxColorValue=255)
    cols[true_unexpressed]<-rgb(44,123,182,maxColorValue=255)
    cols[false_expressed]<-rgb(171,217,233,maxColorValue=255)
    cols[false_unexpressed]<-rgb(253,174,97,maxColorValue=255)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],,col=cols[false_expressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top1000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom1000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top1000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom1000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    write.table(t(c("Top1000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom1000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top1000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom1000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}


#############################################################################################################################################################################
sensitivity_housekeeping<-function(data,output_dir){
    #Check in more detail whether Housekeeping genes are assigned into expressed cluster:
    housekeeping<-read.csv(paste(ref_dir,"/GeneLists/HK_gene_names.txt",sep=""),header=FALSE)


    broad_boundary<-quantile(data$Broad,0.99)
    small_boundary<-quantile(data$Small,0.99)

    housekeeping_in_expressed_cluster<-0
    housekeeping_in_unexpressed_cluster<-0
    data$Housekeeping<-rep("unknown",length(data$Gene))
    for (top in housekeeping$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        data$Housekeeping[index] <- "Housekeeping"
    }
    true_expressed<-which(data$Housekeeping == "Housekeeping" & data$Clustering == "Expressed")
    false_unexpressed<-which(data$Housekeeping == "Housekeeping" & data$Clustering == "Unexpressed")

    cols<-rep(rgb(0.8,0.8,0.8,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(215,25,28,maxColorValue=255)
    cols[false_unexpressed]<-rgb(253,174,97,maxColorValue=255)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Housekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_Housekeeping_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_Housekeeping_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    write.table(t(c("Housekeeping (Eisenberg) in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Houskeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Housekeeping (Eisenberg) in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Houskeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}

#############################################################################################################################################################################
clustering_vs_rma<-function(data,output_dir){
    #check RMA values for genes and plot "Heatmap"
    data$rma<-rep(0,length(data$Gene))
    rma_values<-read.csv(paste(ref_dir,"/NonPregnant_annotated_noChrM_dedup_onlyNM.txt.csv",sep=""),header=TRUE,sep="\t")
    for (gene in rma_values$Gene) {
        index<-which(data$Gene == gene)
        rma_index<-which(rma_values$Gene == gene)
        if (length(index) == 0) {
            next
        }
        if (length(rma_index) == 1) {
            rma_value <- rma_values$Mean[rma_index]
        }
        else  {
            rma_value <- mean(rma_values$Mean[rma_index])
        }
        data$rma[index] <- rep(rma_value,length(index))
    }
    colors<-colorRampPalette(c("green","red"))(3)
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    level<-which(data$rma > 7)
    cols[level]<-colors[3]
    level<-which(data$rma < 4)
    cols[level]<-colors[1]
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Heatmap.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()

    #compare distance to expressed cluster center to RMA
    colors<-colorRampPalette(c("green","red"))(max(floor(data$rma)))
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    for (i in min(floor(data$rma)):max(floor(data$rma))) {
       level<-which(floor(data$rma) == i)
       cols[level]<-colors[i]
    }
}

#############################################################################################################################################################################
clustering_vs_fpkm<-function(data,output_dir) {
    #check FPKM values for genes and plot "Heatmap"
    data$fpkm<-rep(0,length(data$Gene))
    fpkm_values<-read.csv(paste(ref_dir,"/filtered.merged.genes.fpkm_tracking",sep=""),header=TRUE,sep="\t")
    broad_boundary<-quantile(data$Broad,0.99)
    small_boundary<-quantile(data$Small,0.99)
    for (gene in fpkm_values$tracking_id) {
        index<-which(data$Gene == gene)
        fpkm_index<-which(fpkm_values$tracking_id == gene)
        if (length(index) == 0) {
            next
        }
        if (length(fpkm_index) == 1) {
            fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
        }
        else  {
            fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
        }
        data$fpkm[index] <- rep(fpkm_value,length(index))
    }

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Boxplot.png",sep=""))
    boxplot(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
    dev.off()
    true_expressed<-which(data$fpkm > 8 & data$Clustering == "Expressed")
    true_unexpressed<-which(data$fpkm > 0 & data$fpkm < 0.1 & data$Clustering == "Unexpressed")
    false_expressed<-which(data$fpkm > 0 & data$fpkm < 0.1 & data$Clustering ==  "Expressed")
    false_unexpressed<-which(data$fpkm > 8 & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))


    cols<-rep(rgb(0.8,0.8,0.8,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(215,25,28,maxColorValue=255)
    cols[true_unexpressed]<-rgb(44,123,182,maxColorValue=255)
    cols[false_expressed]<-rgb(171,217,233,maxColorValue=255)
    cols[false_unexpressed]<-rgb(253,174,97,maxColorValue=255)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Heatmap.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Highlight.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],col=cols[false_expressed])
    points(data$Broad[true_expressed],data$Small[true_expressed],col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],col=cols[true_unexpressed])
    dev.off()


    write.table(t(c("Top1000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom1000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top1000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom1000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_top100<-function(data,output_dir){
    #Check if snesitivity is higher for Top100 vs Bottom 100:
    top100<-read.csv(paste(ref_dir,"/Top100_NMonly.txt",sep=""),header=FALSE)
    bottom100<-read.csv(paste(ref_dir,"/Bottom100_NMonly.txt",sep=""),header=FALSE)

    broad_boundary<-quantile(data$Broad,0.99)
    small_boundary<-quantile(data$Small,0.99)
    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    for (bottom in bottom100$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom100" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))


    cols<-rep(rgb(0.8,0.8,0.8,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(215,25,28,maxColorValue=255)
    cols[true_unexpressed]<-rgb(44,123,182,maxColorValue=255)
    cols[false_expressed]<-rgb(171,217,233,maxColorValue=255)
    cols[false_unexpressed]<-rgb(253,174,97,maxColorValue=255)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    dev.off()


    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),pch=20)
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],,col=cols[false_expressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom100 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom100 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)

}

##############################################################################################
#use housekeeping and unexpressed genes for evaluation
top1000<-read.csv(paste(ref_dir,"/GeneLists/HK_gene_names.txt",sep=""),header=FALSE)
bottom1000<-read.csv(paste(ref_dir,"/GeneLists/Fantom5_all_lower0.1.txt",sep=""),header=FALSE)

top_in_expressed_cluster<-0
top_in_unexpressed_cluster<-0
bottom_in_expressed_cluster<-0
bottom_in_unexpressed_cluster<-0
data$Expressed<-rep("unknown",length(data$Gene))
for (top in top1000$V1) {
    index<-which(data$Gene == top)
    if (length(index) == 0) {
        next
    }
    data$Expressed[index] <- "Top1000"
}
for (bottom in bottom1000$V1) {
    index<-which(data$Gene == bottom)
    if (length(index) == 0) {
        next
    }
    data$Expressed[index] <- "Bottom1000"
}

data$Prediction_Expr<-rep(0,length(data$Gene))
data$Prediction_Unexpr<-rep(0,length(data$Gene))
for (i in 1:1000) {
    set.seed(i)
    data$Training<-rep(FALSE,length(data$Gene))
    training_top_index<-sample(which(data$Expressed == "Top1000"),300)
    training_bottom_index<-sample(which(data$Expressed == "Bottom1000"),300)
    data$Training[training_top_index]<-TRUE
    data$Training[training_bottom_index]<-TRUE
    training_data<-data.frame(Broad=c(data$Broad[training_top_index],data$Broad[training_bottom_index]),Small=c(data$Small[training_top_index],data$Small[training_bottom_index]))
    type<-rep("Unexpressed",length(training_data$Broad))
    type[1:300]<-"Expressed"
    training_data$Type<-as.factor(type)
    svm_model<-svm(Type ~ .,data=training_data,type="C-classification")

    test_data_top_index<-which(data$Training == FALSE)
    test_data_bottom_index<-which(data$Training == FALSE)
    test_data<-data.frame(Broad=c(data$Broad[test_data_top_index],data$Broad[test_data_bottom_index]),Small=c(data$Small[test_data_top_index],data$Small[test_data_bottom_index]),Index=c(test_data_top_index,test_data_bottom_index))
    type<-rep("Unexpressed",length(test_data$Broad))
    type[1:length(test_data_top_index)]<-"Expressed"
    test_data$Type<-as.factor(type)
    prediction<-predict(svm_model,test_data)
    data$Prediction_Expr[test_data$Index[which(prediction == "Expressed")]]<-data$Prediction_Expr[test_data$Index[which(prediction == "Expressed")]]+1
    data$Prediction_Unexpr[test_data$Index[which(prediction == "Unexpressed")]]<-data$Prediction_Unexpr[test_data$Index[which(prediction == "Unexpressed")]]+1
}

true_expressed<-which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Top1000")
false_unexpressed<-which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Top1000")
true_unexpressed<-which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Bottom1000")
false_expressed<-which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Bottom1000")


col=rep(rgb(0,0,0,0.2),length(data$Gene))
col[true_expressed]<-"red"
col[true_unexpressed]<-"green"
col[false_unexpressed]<-"blue"
col[false_expressed]<-"blue"
png(paste(output_dir,"/Prediction_Top1000.png",sep=""))
broad_boundary<-quantile(data$Broad,0.99)
small_boundary<-quantile(data$Small,0.99)
plot(data$Broad,data$Small,col=col,xlim=c(0,broad_boundary),ylim=c(0,small_boundary),xlab="Broad Coverage",ylab="Small Coverage")
dev.off()

length(which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Top1000"))
length(which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Top1000"))
length(which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Bottom1000"))
length(which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Bottom1000"))

data$Clustering<-rep("Unexpressed",length(data$Gene))
expressed_index<-which(data$Prediction_Expr > data$Prediction_Unexpr)
data$Clustering[expressed_index]<-"Expressed"


#############################################################################################
# Only analyze genes which are consistently called in >75% of the SVM runs
expr_consistent<-which(data$Prediction_Expr > 0.75 * (data$Prediction_Expr+data$Prediction_Unexpr)) 
unexpr_consistent<-which(data$Prediction_Unexpr > 0.75 * (data$Prediction_Expr+data$Prediction_Unexpr)) 
data$Clustering<-rep("unknown",length(data$Gene))
data$Clustering[expr_consistent]<-"Expressed"
data$Clustering[unexpr_consistent]<-"Unexpressed"

col=rep(rgb(0,0,0,0.2),length(data$Gene))
col[expr_consistent]<-"red"
col[unexpr_consistent]<-"green"
png(paste(output_dir,"/Consistent_calls.png",sep=""))
broad_boundary<-quantile(data$Broad,0.99)
small_boundary<-quantile(data$Small,0.99)
plot(data$Broad,data$Small,col=col,xlim=c(0,broad_boundary),ylim=c(0,small_boundary))
dev.off()

col=rep("blue",length(data$Gene))
col[expr_consistent]<-rgb(0,0,0,0.2)
col[unexpr_consistent]<-rgb(0,0,0,0.2)
png(paste(output_dir,"/Inconsistent_calls.png",sep=""))
plot(data$Broad,data$Small,col=col,xlim=c(0,broad_boundary),ylim=c(0,small_boundary))
dev.off()

analyze_and_plot_clustering(data,output_dir)
sensitivity_top5000(data,output_dir)
sensitivity_top1000(data,output_dir)
sensitivity_top100(data,output_dir)
clustering_vs_rma(data,output_dir)
clustering_vs_fpkm(data,output_dir)

# Check distribution of Top1000 and Bottom1000 Genes
top_bottom<-which((data$Expressed == "Top1000" | data$Expressed == "Bottom1000") & (data$Broad < 25) & (data$Small < 20))
top<-which(data$Expressed[top_bottom] == "Top1000")
bottom<-which(data$Expressed[top_bottom] == "Bottom1000")
col=rep("black",length(data$Broad[top_bottom]))
col[top]<-"red"
col[bottom]<-"green"
png(paste(output_dir,"/distribution_top_and_bottom1000.png",sep=""))
plot(data$Broad[top_bottom],data$Small[top_bottom],col=col,xlim=c(0,broad_boundary),ylim=c(0,small_boundary))
dev.off()
#Kernel density estimation
library(MASS)
kd_estimation<-kde2d(data$Broad[top_bottom],data$Small[top_bottom],n=50)

png(paste(output_dir,"/KDE_Contour.png",sep=""))
contour(kd_estimation,xlab="Broad Normalized Coverage",ylab="Small Coverage",xlim=c(0,broad_boundary),ylim=c(0,small_boundary))
dev.off()

png(paste(output_dir,"/Heatmap.png",sep="")) 
image(kd_estimation,xlab="Broad Normalized Coverage",ylab="Small Coverage",xlim=c(0,broad_boundary),ylim=c(0,small_boundary))
dev.off()
png(paste(output_dir,"/Perspective.png",sep="")) 
persp(kd_estimation, phi = 45, theta = 30,xlab="Broad Normalized Coverage",ylab="Small Coverage",zlab="Density",xlim=c(0,broad_boundary),ylim=c(0,small_boundary))
dev.off()

png(paste(output_dir,"/Heatmap_contour.png",sep=""))
# fancy contour with image
image(kd_estimation,xlab="Broad Normalized Coverage",ylab="Small Coverage",xlim=c(0,broad_boundary),ylim=c(0,small_boundary)); contour(kd_estimation, add = T)
dev.off()

# fancy perspective
png(paste(output_dir,"/Perspective_angle.png",sep=""))
persp(kd_estimation, phi = 30, theta = 35, shade = .1, border = NA,xlab="Broad Normalized Coverage",ylab="Small Coverage",zlab="Density")
dev.off()



#Check up whether expressed genes have higher FPKM
data$fpkm<-rep(0,length(data$Gene))
fpkm_values<-read.csv(paste(ref_dir,"/filtered.merged.genes.fpkm_tracking",sep=""),header=TRUE,sep="\t")
for (gene in fpkm_values$tracking_id) {
    index<-which(data$Gene == gene)
    fpkm_index<-which(fpkm_values$tracking_id == gene)
    if (length(index) == 0) {
        next
    }
    if (length(fpkm_index) == 1) {
        fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
    }
    else  {
        fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
    }
        data$fpkm[index] <- rep(fpkm_value,length(index))
}
w_cons<-wilcox.test(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
t_cons<-t.test(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
med_expressed_cons<-median(data$fpkm[which(data$Clustering=="Expressed")])
med_unexpressed_cons<-median(data$fpkm[which(data$Clustering=="Unxpressed")])
sd_expressed_cons<-sd(data$fpkm[which(data$Clustering=="Expressed")])
sd_unexpressed_cons<-sd(data$fpkm[which(data$Clustering=="Unxpressed")])

png(paste(output_dir,"Broad_vs_Small_TSS_coverage_FPKM_Boxplot.png",sep=""))
boxplot(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")],ylim=c(0,100),names=c("Expressed","Unexpressed"),col=c("red","green"))
dev.off()

#Check up whether predicted expressed genes have higher RMA
data$rma<-rep(0,length(data$Gene))
rma_values<-read.csv(paste(ref_dir,"/NonPregnant_annotated_noChrM_dedup_onlyNM.txt.csv",sep=""),header=TRUE,sep="\t")
for (gene in rma_values$Gene) {
    index<-which(data$Gene == gene)
    rma_index<-which(rma_values$Gene == gene)
    if (length(index) == 0) {
        next
    }
    if (length(rma_index) == 1) {
        rma_value <- rma_values$Mean[rma_index]
    }
    else  {
        rma_value <- mean(rma_values$Mean[rma_index])
    }
    data$rma[index] <- rep(rma_value,length(index))
}
w_cons_rma<-wilcox.test(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")])
t_cons_rma<-t.test(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")])
med_expressed_cons_rma<-median(data$rma[which(data$Clustering=="Expressed")])
med_unexpressed_cons_rma<-median(data$rma[which(data$Clustering=="Unexpressed")])

png(paste(output_dir,"/Broad_vs_Small_TSS_coverage_RMA_Boxplot.png",sep=""))
boxplot(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")],names=c("Expressed","Unexpressed"),col=c("red","green"))
dev.off()

