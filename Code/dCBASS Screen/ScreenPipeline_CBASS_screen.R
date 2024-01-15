library(graphics)
library(stats)
library(Hmisc)
library(zoo)
library(gplots)
library(corrplot)
library(gtools)
library(caTools)
library(chron)

#====== Set directories, colors, shapes and cut-off ========
here_path = file.choose()
Load_dir = paste0(here_path,"Plates/") 
Out_dir = paste0(here_path,"Out/")
clr = c("#0B6199","#EF4D04","#03A266","#EF8704")
cols_main = c("#0B6199","#EF4D04","#03A266","#EF8704")
cols_light = c("#4C8BB5","#FF9361","#48BC90","#FFB961")

rep_pch = c(15,19,17)
rep_pch_out = c(0,1,2)

#Define cutoff for Residuals to call hits according to mean of wt vs wt baseline comparison
cutoff = 0.2

#====== Set functions ========
#Like grep, but returns only exact match
grep_exact <- function(x,v)
{
  match_vec = grep(x, v, fixed = TRUE)
  uniq_match = c()
  
  #To make sure that well_index points to the exact one and only one well
  if (length(match_vec) > 0)
  {
    for(i in 1:length(match_vec))
    {
      match_vec_temp = match_vec[i]
      if (nchar(match_vec_temp) > 0 && x == v[match_vec_temp])
      {uniq_match = c(uniq_match, match_vec_temp)}
    }
  }
  return(uniq_match)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="complete"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#====== Vector files to re-formated ========

#Plates
All_files = list.files(Load_dir)
Plates = All_files[grep("txt",All_files)]
Plates = substr(Plates,1,nchar(Plates)-4) #gets rid of the ".txt"

#PlateDatabase
file_id = paste0(here_path,"R-scripts/Database.txt",collapse=NULL)
PlateDatabase = read.table(file_id,header = T, sep = "\t")
file_id = paste0(here_path,"R-scripts/Map.txt",collapse=NULL)
Map = read.table(file_id,header = T, sep = "\t")

rm(file_id)

#====== Read data from plates into lists ======
OD_data = list()
Plate_names = c()

for(p in 1:length(Plates))
{
  file_id = paste0(Load_dir,Plates[p],".txt")
  plate_data = read.table(file = file_id,sep = "\t",header = T)
  
  detection = PlateDatabase$Detection[match(Plates[p],substr(PlateDatabase$file_id,1,nchar(PlateDatabase$file_id)-4))]
  plate_name = substr(Plates[p],1,nchar(Plates[p])-3)
  Plate_names = c(Plate_names,plate_name)
  OD_data[[length(OD_data)+1]] = plate_data
}
#rename plates to match the Plate_nr in Plate database
Plates = Plate_names
names(OD_data) = Plate_names

#====== Plot automatically for OD of all plates ===============
Out_dir = paste0(here_path,"Out/")
uniqdrugs= as.character(unique(Map$Drug))

file_id = paste0(Out_dir,"OD_raw_all_plates.pdf",collapse=NULL)
pdf(file_id, useDingbats = T)
for(p in 1:length(OD_data))
{
  
  OD_plate = OD_data[[p]] 
  Plate_id = names(OD_data)[p]
  
  par(mfrow=c(8,12),mar=c(1,0.5,1,0.5),cex=0.32)
  
  uniqdrugs= as.character(unique(Map$Drug))
  
  for(drug in uniqdrugs)
  {
    drug_wells = as.character(Map$Well[grep_exact(drug, Map$Drug)])
    
    drug_ConcMock = as.character(Map$ConcMock[match(drug_wells, Map$Well)])
    
    drug_data = OD_plate[match(drug_wells ,names(OD_plate))]
    
    plot(OD_plate[[1]],drug_data[[1]], ylim= c(0,1), xlim = c(0,12), col=clr[1], pch=19, main=drug,axes=F,frame=T)
    for(i in 2:length(drug_data)) 
    {points(OD_plate[[1]],drug_data[[i]], col= clr[i], pch=19)}
    legend("topleft", legend = paste(drug_wells, drug_ConcMock), col = clr, pch=19, bty="n")
    
  }
  text(x = 9,y = 0.1,labels = Plate_id,cex = 2)
  
}
dev.off()

#====== Background correction of OD for all plates plus plots ======

file_id = paste0(Out_dir,"OD_corrected_all_plates.pdf",collapse=NULL)
pdf(file_id, useDingbats = T)

OD_data_corrected = list()

for(p in 1:length(OD_data))
{
  
  OD_plate = OD_data[[p]] 
  Plate_id = names(OD_data)[p]
  
  bg_data= OD_plate[1,]
  bg_mean= t(as.data.frame(sapply(bg_data,FUN=mean)))
  OD_plate_corrected = as.data.frame(t(as.data.frame(t(OD_plate)) - bg_mean))
  OD_plate_corrected[1] = OD_plate[1]
  
  OD_data_corrected[[length(OD_data_corrected)+1]] = OD_plate_corrected
  
  par(mfrow=c(8,12),mar=c(1,0.5,1,0.5),cex=0.32)
  
  uniqdrugs= as.character(unique(Map$Drug))
  
  for(drug in uniqdrugs)
  {
    drug_wells = as.character(Map$Well[grep_exact(drug, Map$Drug)])
    
    drug_ConcMock = as.character(Map$ConcMock[match(drug_wells, Map$Well)])
    
    drug_data = OD_plate_corrected[match(drug_wells ,names(OD_plate_corrected))]
    
    plot(OD_plate_corrected[[1]],drug_data[[1]], ylim= c(0,1), col=clr[1], pch=19, main=drug,axes=F,frame=T)
    for(i in 2:length(drug_data)) 
    {points(OD_plate_corrected[[1]],drug_data[[i]], col= clr[i], pch=19)}
    legend("topleft", legend = paste(drug_wells, drug_ConcMock), col = clr, pch=19, bty="n")
    
  }
  text(x = 9,y = 0.1,labels = Plate_id,cex = 2)
  
}
dev.off()
names(OD_data_corrected) = Plates

#====== Calculate AUC of corrected ODs for all plates ========

for(p in 1:length(OD_data_corrected))
{
  
  OD_plate_corrected = OD_data_corrected[[p]]
  Plate_id = names(OD_data_corrected)[p]
  
  OD_AUC = sapply(OD_plate_corrected, FUN= trapz, x = OD_plate_corrected[[1]])
  OD_AUC = as.data.frame(OD_AUC[-c(1)])
  
  names(OD_AUC) = Plate_id
  
  if(p==1)
  {All_OD_AUC = OD_AUC} else
  {All_OD_AUC = as.data.frame(cbind(All_OD_AUC, OD_AUC))}
  
}


#====== Open plotting file ==============
file_id = paste0(Out_dir,"SummaryScreen_forIllust.pdf",collapse=NULL)
pdf(file_id, useDingbats = F)
par(family = "sans")

#====== Center data to ensure safer comparability ==============
boxplot(All_OD_AUC,pch=rep(rep_pch,each=3),border = cols_light[1],col="white",main="Raw AUC")
median_All_OD_AUC = sapply(All_OD_AUC,FUN=median,na.rm=T)

centered_All_OD_AUC = as.data.frame(t(as.data.frame(t(All_OD_AUC))/median_All_OD_AUC))
boxplot(centered_All_OD_AUC,pch=rep(rep_pch,each=3),border = cols_light[1],col="white",main="Centered AUC")

All_OD_AUC = centered_All_OD_AUC

#====== Show high replicate correlation ==============
pairs(All_OD_AUC[c(1:3,4:6)],lower.panel = panel.cor,pch=19,
      col=alpha(cols_light[1],alpha=0.5),
      main="Replicate correlation")

#====== Baseline analysis - WT =======

plot(0,0,type="n",axes=F,frame=F,ylab="",xlab="")
text(-0.1,0,labels = "Baseline (WT) analysis")

#Plot replicate correlation to get a baseline residual distance 
r_data = All_OD_AUC[grep("wt",names(All_OD_AUC))] #WT data
par(mfrow=c(2,2),mar=c(4,4,4,4))
for(r in 2:3)
{
  lm_rep_r = lm(r_data[[1]]~r_data[[r]]-1)
  res_r = as.data.frame(lm_rep_r$residuals)
  row.names(res_r) = row.names(r_data)
  names(res_r) = paste0("WTrep 1x",r)
  Rsq = 1-(sum(res_r^2)/sum((r_data[[1]]-mean(r_data[[r]],na.rm=T))^2))
  
  highlight_col = "grey20"
  
  wells = row.names(res_r)[grep(T,res_r[[1]]>cutoff)]
  drugs = Map$Drug[match(wells,Map$Well)]
  
  plot(r_data[[r]],r_data[[1]],main=paste0("rep ",r),ylab="AUC WT replicate i",xlab="AUC WT replicate j",
       pch=rep_pch[r],col="gray",ylim=c(0,1.5))
  #points(r_data[[r]],r_data[[1]],pch=1,col="gray50")
  abline(lm_rep_r,col="gray50")
  points(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)],
         pch=rep_pch[r],col=alpha(highlight_col,alpha=0.5))
  text(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)]+0.05,
       labels=drugs,cex=0.5,col=highlight_col)
  legend("topleft",legend = c(paste0("slope=", round(lm_rep_r$coefficients[[1]],3)),
                              paste0("Rsq=", round(Rsq[[1]],3))),bty="n")
  
  if(r==2)
  {Resids_baseline = res_r} else
  {Resids_baseline = as.data.frame(cbind(Resids_baseline,res_r))}
  
  if(r==2)
  {
    Rsqs=Rsq
    Slopes = round(lm_rep_r$coefficients[[1]],3)
  } else
  {
    Rsqs=c(Rsqs,Rsq)
    Slopes = c(Slopes,round(lm_rep_r$coefficients[[1]],3))
  }
  
  if(r==3) # get correlation between rep 2 & 3
  {
    lm_rep_r = lm(r_data[[2]]~r_data[[r]]-1)
    res_r = as.data.frame(lm_rep_r$residuals)
    row.names(res_r) = row.names(r_data)
    names(res_r) = paste0("WTrep 2x",r)
    Rsq = 1-(sum(res_r^2)/sum((r_data[[2]]-mean(r_data[[r]],na.rm=T))^2))
    
    wells = row.names(res_r)[grep(T,res_r[[1]]>cutoff)]
    drugs = Map$Drug[match(wells,Map$Well)]
    
    plot(r_data[[r]],r_data[[2]],main=paste0("rep ",r),ylab="AUC WT replicate i",xlab="AUC WT replicate j",
         pch=rep_pch[r],col="gray")
    #points(r_data[[r]],r_data[[2]],pch=1,col="gray50")
    abline(lm_rep_r,col="gray50")
    points(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[2]][grep(T,res_r[[1]]>cutoff)],
           pch=rep_pch[r],col=alpha(highlight_col,alpha=0.7))
    text(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[2]][grep(T,res_r[[1]]>cutoff)]+0.05,
         labels=drugs,cex=0.6,col=highlight_col)
    legend("topleft",legend = c(paste0("slope=", round(lm_rep_r$coefficients[[1]],3)),
                                paste0("Rsq=", round(Rsq[[1]],3))),bty="n")
    
    Resids_baseline = as.data.frame(cbind(Resids_baseline,res_r))
    
    
    Rsqs=c(Rsqs,Rsq)
    Slopes = c(Slopes,round(lm_rep_r$coefficients[[1]],3))
  }
}
boxplot(Resids_baseline,col ="white", border = "grey50",pch=rep_pch,main="Residuals baseline")
points(c(-1,10),c(cutoff,cutoff),type="l")

#overlay replicates
par(mfrow=c(1,1),mar=c(4,4,4,4))
r_data = All_OD_AUC[grep("wt",names(All_OD_AUC))] #WT data
for(r in 2:3)
{
  lm_rep_r = lm(r_data[[1]]~r_data[[r]]-1)
  res_r = as.data.frame(lm_rep_r$residuals)
  row.names(res_r) = row.names(r_data)
  names(res_r) = paste0("WTrep 1x",r)
  Rsq = 1-(sum(res_r^2)/sum((r_data[[1]]-mean(r_data[[r]],na.rm=T))^2))
  
  highlight_col = "grey50"
  
  wells = row.names(res_r)[grep(T,res_r[[1]]>cutoff)]
  drugs = Map$Drug[match(wells,Map$Well)]
  
  if(r==2)
  {
    plot(r_data[[r]],r_data[[1]],main="WT residuals baseline",ylab="AUC replicate i",xlab="AUC replicate j",
         pch=rep_pch[r],col=alpha(highlight_col,alpha=0.2),ylim=c(0,1.5))
  } else
  {points(r_data[[r]],r_data[[1]],pch=rep_pch[r],col=alpha(highlight_col,alpha=0.2))}
  
  #points(r_data[[r]],r_data[[1]],pch=1,col=alpha(highlight_col,alpha=0.7))
  abline(lm_rep_r,col="gray50")
  points(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)],
         pch=rep_pch[r],col=alpha(highlight_col,alpha=1))
  text(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)]+0.05,
       labels=drugs,cex=0.6,col=highlight_col)
  
  if(r==3) # get correlation between rep 2 & 3
  {
    lm_rep_r = lm(r_data[[2]]~r_data[[r]]-1)
    res_r = as.data.frame(lm_rep_r$residuals)
    row.names(res_r) = row.names(r_data)
    names(res_r) = paste0("WTrep 2x",r)
    Rsq = 1-(sum(res_r^2)/sum((r_data[[2]]-mean(r_data[[r]],na.rm=T))^2))
    
    wells = row.names(res_r)[grep(T,res_r[[1]]>cutoff)]
    drugs = Map$Drug[match(wells,Map$Well)]
    
    points(r_data[[r]],r_data[[2]],pch=rep_pch[r],col=alpha(highlight_col,alpha=0.2))
    #points(r_data[[r]],r_data[[2]],pch=1,col=alpha(highlight_col,alpha=0.7))
    abline(lm_rep_r,col="gray50")
    points(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[2]][grep(T,res_r[[1]]>cutoff)],
           pch=rep_pch[r],col=alpha(highlight_col,alpha=1))
    text(r_data[[r]][grep(T,res_r[[1]]>cutoff)],r_data[[2]][grep(T,res_r[[1]]>cutoff)]+0.05,
         labels=drugs,cex=0.6,col=highlight_col)
  }
}
legend("topleft", legend = paste0("Slope=",Slopes),col=highlight_col,pch=rep_pch,
       bty="n")
legend("bottomright", legend = paste0("Rsq=",round(Rsqs,3)),col=highlight_col,pch=rep_pch,
       bty="n")

boxplot(Resids_baseline,col ="white", border = "grey",pch=19,main="WT residuals baseline",ylab="Residuals")
points(c(-1,10),c(cutoff,cutoff),type="l")

#====== Calculate mean residuals for baseline ======
mean_resids = as.data.frame(sapply(as.data.frame(t(Resids_baseline)),FUN=mean,na.rm=T))
names(mean_resids) = "mean_resids"
sd_resids = as.data.frame(sapply(as.data.frame(t(Resids_baseline)),FUN=sd,na.rm=T))
names(sd_resids) = "sd_resids"
Resids_baseline = as.data.frame(cbind(Resids_baseline,mean_resids,sd_resids))
new_names = paste0(Map$Drug[match(row.names(Resids_baseline),Map$Well)],
                   Map$ConcMock[match(row.names(Resids_baseline),Map$Well)])
row.names(Resids_baseline) = new_names

Sorted_resids_baseline = Resids_baseline[rev(order(Resids_baseline$mean_resids)),]
barplot(Sorted_resids_baseline$mean_resids,ylab="Residuals baseline",xlab="Drug",col = "grey",border = "grey52",
        ylim=c(-0.1,0.5))
plot(density(Resids_baseline$mean_resids),lwd=2,col="grey52")
points(c(-1,400),c(cutoff,cutoff),type="l")

#Zoom in on the top 30 drugs
bar_cols = c(rep("grey",length(Sorted_resids_baseline[[1]])))
bar_cols[grep(T,Sorted_resids_baseline$mean_resids>=cutoff)] = alpha("gray52",alpha=1)

barplot(c(Sorted_resids_baseline$mean_resids[1:30],mean(Resids_baseline$mean_resids)),
        names.arg = c(row.names(Sorted_resids_baseline)[1:30],"mean of residuals"),
        ylab="Residuals",xlab="Drug",border = "grey52",las=2,cex.names = 0.5,
        col=bar_cols,ylim=c(0,0.5),main= "Residuals baseline")
points(c(-1,100),c(cutoff,cutoff),type="l")

#====== CBASS vs WT Analysis =======

plot(0,0,type="n",axes=F,frame=F,ylab="",xlab="")
text(-0.1,0,labels = "CBASS analysis")

#Plot linear fits & call hits based on Residual cutoff 
#plot data as residuals
par(mfrow=c(2,2),mar=c(4,4,4,4))
for(r in 1:3)
{
  r_data = All_OD_AUC[grep(r,names(All_OD_AUC))]
  lm_rep_r = lm(r_data[[1]]~r_data[[2]]-1)
  res_r = as.data.frame(lm_rep_r$residuals)
  row.names(res_r) = row.names(r_data)
  names(res_r) = paste0("rep ",r)
  Rsq = 1-(sum(res_r^2)/sum((r_data[[2]]-mean(r_data[[2]],na.rm=T))^2))
  highlight_col = cols_main[r]
  wells = row.names(res_r)[grep(T,res_r[[1]]>cutoff)]
  drugs = Map$Drug[match(wells,Map$Well)]
  
  plot(r_data[[2]],r_data[[1]],main=paste0("rep ",r),ylab="AUC deltaCBASS",xlab="AUC WT",
       pch=rep_pch[r],col="gray",ylim=c(0,1.5),xlim=c(0,1.5))
  abline(lm_rep_r,col="gray50")
  points(r_data[[2]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)],
         pch=rep_pch[r],col=alpha(highlight_col,alpha=1))
  text(r_data[[2]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)]+0.05,
       labels=drugs,cex=0.5,col=highlight_col)
  legend("topleft",legend = c(paste0("slope=", round(lm_rep_r$coefficients[[1]],3)),
                              paste0("Rsq=", round(Rsq[[1]],3))),bty="n")
  if(r==1)
  {Resids = res_r} else
  {Resids = as.data.frame(cbind(Resids,res_r))}
  
  if(r==1)
  {
    Rsqs=Rsq
    Slopes = round(lm_rep_r$coefficients[[1]],3)
  } else
  {
    Rsqs=c(Rsqs,Rsq)
    Slopes = c(Slopes,round(lm_rep_r$coefficients[[1]],3))
  }
  
}
boxplot(Resids,col ="white", border = alpha(cols_light,alpha=0.5),pch=19,main="Residuals")
points(c(-1,10),c(cutoff,cutoff),type="l")

#overlay replicates
par(mfrow=c(1,1),mar=c(4,4,4,4))
for(r in 1:3)
{
  r_data = All_OD_AUC[grep(r,names(All_OD_AUC))]
  lm_rep_r = lm(r_data[[1]]~r_data[[2]]-1)
  res_r = as.data.frame(lm_rep_r$residuals)
  row.names(res_r) = row.names(r_data)
  names(res_r) = paste0("rep ",r)
  
  highlight_col = cols_light[r]
  
  wells = row.names(res_r)[grep(T,res_r[[1]]>cutoff)]
  drugs = Map$Drug[match(wells,Map$Well)]
  
  if(r==1)
  {
    plot(r_data[[2]],r_data[[1]],main="CBASS vs WT residuals",ylab="AUC deltaCBASS",xlab="AUC WT",
         pch=rep_pch[r],col="gray70",ylim=c(0,1.4),xlim=c(0,1.4))
  } else
  {points(r_data[[2]],r_data[[1]],pch=rep_pch[r],col="gray70")}
  abline(lm_rep_r,col="gray50")
  abline(a = 0.2, b = 1)
  points(r_data[[2]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)],
         pch=19,col=alpha(highlight_col,alpha=0.7))
  points(r_data[[2]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)],
         pch=rep_pch[r],col=alpha(highlight_col,alpha=1))
  
  text(r_data[[2]][grep(T,res_r[[1]]>cutoff)],r_data[[1]][grep(T,res_r[[1]]>cutoff)]+0.05,
       labels=drugs,cex=0.6,col=highlight_col)
  
}
legend("topleft", legend = paste0("Slope=",Slopes),col=cols_light,pch=19,
       bty="n")
legend("bottomright", legend = paste0("Rsq=",round(Rsqs,3)),col=cols_light,pch=19,
       bty="n")
boxplot(Resids,col ="white", border = alpha(cols_light,alpha=0.5),pch=19,main="Residuals")
points(c(-1,10),c(cutoff,cutoff),type="l")

#====== Calculate mean residuals CBASS ======

mean_resids = as.data.frame(sapply(as.data.frame(t(Resids)),FUN=mean,na.rm=T))
names(mean_resids) = "mean_resids"
sd_resids = as.data.frame(sapply(as.data.frame(t(Resids)),FUN=sd,na.rm=T))
names(sd_resids) = "sd_resids"
Resids = as.data.frame(cbind(Resids,mean_resids,sd_resids))
new_names = paste0(Map$Drug[match(row.names(Resids),Map$Well)],
                   Map$ConcMock[match(row.names(Resids),Map$Well)])
row.names(Resids) = new_names

#====== Get actual residuals vs baseline residuals =========

plot(Resids$mean_resids,Resids_baseline$mean_resids,col=alpha("grey",alpha=1),
     pch=19,main="Cutoff comparison",ylim=c(-0.1,0.5),xlim=c(-0.1,0.5),cex=1.5)
points(c(cutoff,cutoff),c(-1,10),type="l")
points(c(-1,10),c(cutoff,cutoff),type="l")
text(Resids$mean_resids[grep(T,Resids$mean_resids>cutoff)],
     Resids_baseline$mean_resids[grep(T,Resids$mean_resids>cutoff)] + 0.01,
     labels=row.names(Resids)[grep(T,Resids$mean_resids>cutoff)],cex=0.5)

#====== Get sorted residuals ========

Sorted_resids = Resids[rev(order(Resids$mean_resids)),]
barplot(Sorted_resids$mean_resids,ylab="Residuals",xlab="Drug",col = "grey",border = "grey52",
        ylim=c(-0.1,0.5))
plot(density(Resids$mean_resids),lwd=2,col="grey52")

#Zoom in on the top 30 drugs
bar_cols = c(rep("grey",length(Sorted_resids[[1]])))
bar_cols[grep(T,Sorted_resids$mean_resids>=cutoff)] = alpha(cols_light[1],alpha=0.7)

barplot(c(Sorted_resids$mean_resids[1:30],mean(Resids$mean_resids)),
        names.arg = c(row.names(Sorted_resids)[1:30],"mean of residuals"),
        ylab="Residuals",xlab="Drug",border = "grey52",las=2,cex.names = 0.5,
        col=bar_cols,ylim=c(0,0.5))
points(c(-1,100),c(cutoff,cutoff),type="l")

Sorted_resids = Sorted_resids[grep(T,Sorted_resids$mean_resids>=cutoff),]

#====== Calculate p-values =========
for(w in 1:length(Sorted_resids[[1]]))
{
  w_res = as.numeric(as.matrix(Sorted_resids[w,1:3]))
  w_res_baseline = Resids_baseline[match(row.names(Sorted_resids)[w],row.names(Resids_baseline)),]
  w_res_baseline = abs(as.numeric(as.matrix(w_res_baseline[1:3])))
  test = t.test(w_res,w_res_baseline)
  w_p_val = round(test$p.value,3)
  
  if(w==1)
  {P_vals = c()}
  P_vals = c(P_vals,w_p_val)
}

for(w in 1:length(Resids[[1]]))
{
  w_res_all = as.numeric(as.matrix(Resids[w,1:3]))
  w_res_baseline_all = Resids_baseline[match(row.names(Resids)[w],row.names(Resids_baseline)),]
  w_res_baseline_all = abs(as.numeric(as.matrix(w_res_baseline_all[1:3])))
  test_all = t.test(w_res_all,w_res_baseline_all)
  w_p_val_all = round(test_all$p.value,3)
  if(w==1)
  {P_vals_all = c()}
  P_vals_all = c(P_vals_all,w_p_val_all)
}
Sorted_resids = as.data.frame(cbind(Sorted_resids,P_vals))
names(Sorted_resids)[(length(Sorted_resids))] = "p_vals"
P_vals = p.adjust(P_vals,method = "BH")
Sorted_resids = as.data.frame(cbind(Sorted_resids,P_vals))

#get mean(abs(residuals WT))
dummy = abs(Resids_baseline[match(row.names(Sorted_resids),
                                  row.names(Resids_baseline)),])
dummy[4] = as.data.frame(sapply(as.data.frame(t(dummy[1:3])),FUN=mean,na.rm=T))
dummy[5] = as.data.frame(sapply(as.data.frame(t(dummy[1:3])),FUN=sd,na.rm=T))

Sorted_resids = as.data.frame(cbind(Sorted_resids,dummy[4:5]))
names(Sorted_resids)[8:9] = paste0(names(Sorted_resids)[8:9],"_baseline") 

#====== Get ggplots for barplot with error bars =========

Sorted_resids_ggplots = as.data.frame(Sorted_resids$mean_resids)
Sorted_resids_ggplots = as.data.frame(cbind(Sorted_resids_ggplots,
                                            Sorted_resids$sd_resids))
Sorted_resids_ggplots = as.data.frame(cbind(Sorted_resids_ggplots,
                                            Sorted_resids$P_vals))
Sorted_resids_ggplots = as.data.frame(cbind(Sorted_resids_ggplots,
                                            c(1:length(Sorted_resids[[1]]))))
Sorted_resids_ggplots = as.data.frame(cbind(Sorted_resids_ggplots,
                                            row.names(Sorted_resids)))
Sorted_resids_ggplots = as.data.frame(cbind(Sorted_resids_ggplots,
                                            rep("CBASS",length(Sorted_resids[[1]]))))

names(Sorted_resids_ggplots) = c("mean_resids","sd_resid","P_vals","name","drug","background")

#add the baseline to Sorted_resids_ggplots

Sorted_resids_ggplots_baseline = as.data.frame(Sorted_resids$mean_resids_baseline)
Sorted_resids_ggplots_baseline = as.data.frame(cbind(Sorted_resids_ggplots_baseline,
                                                     Sorted_resids$sd_resids_baseline))
Sorted_resids_ggplots_baseline = as.data.frame(cbind(Sorted_resids_ggplots_baseline,
                                                     Sorted_resids$P_vals))
Sorted_resids_ggplots_baseline = as.data.frame(cbind(Sorted_resids_ggplots_baseline,
                                                     c(1:length(Sorted_resids[[1]]))))
Sorted_resids_ggplots_baseline = as.data.frame(cbind(Sorted_resids_ggplots_baseline,row.names(Sorted_resids)))
Sorted_resids_ggplots_baseline = as.data.frame(cbind(Sorted_resids_ggplots_baseline,
                                                     rep("WT",length(Sorted_resids[[1]]))))

names(Sorted_resids_ggplots_baseline) = c("mean_resids","sd_resid","P_vals","name","drug","background")

#Fuse orted_resids_ggplots
Sorted_resids_ggplots = as.data.frame(rbind(Sorted_resids_ggplots,Sorted_resids_ggplots_baseline))
rm(Sorted_resids_ggplots_baseline,dummy)

p3 <- ggplot(data = Sorted_resids_ggplots, aes(x=drug, y=mean_resids, fill=background)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_resids-sd_resid, ymax=mean_resids+sd_resid), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=c(cols_light[1],'grey60'))

p3 + ggtitle(paste0("Screen hits, cutoff =",cutoff)) + ylim(0, 0.5)

#Output p_vals
plot(0,0,col="white",main="Output p_vals")
legend("topleft",legend = paste0(row.names(Sorted_resids),"_",c(1:length(Sorted_resids[1])),
                                 " p_val =",Sorted_resids$P_vals))
legend("bottomright",legend = paste0(row.names(Sorted_resids),"_",c(1:length(Sorted_resids[1])),
                                     " p_val =",Sorted_resids$p_vals))

#====== Get ggplots for with points to overlay in barplot =========

Sorted_resis_points = as.data.frame(cbind(c(Sorted_resids$`rep 1`,Sorted_resids$`rep 2`,Sorted_resids$`rep 3`),
                                          rep(row.names(Sorted_resids),3),
                                          rep("CBASS",length(Sorted_resids[[1]])*3),
                                          rep(seq(1,length(Sorted_resids[[1]])),3)))

baseline_dummy = Sorted_resids_baseline[match(row.names(Sorted_resids),row.names(Sorted_resids_baseline)),]


Sorted_resis_points = as.data.frame(rbind(Sorted_resis_points,
                                          as.data.frame(cbind(abs(c(baseline_dummy$`WTrep 1x2`,baseline_dummy$`WTrep 1x3`,baseline_dummy$`WTrep 2x3`)),
                                                              rep(row.names(baseline_dummy),3),
                                                              rep("WT",length(baseline_dummy[[1]])*3),
                                                              rep(seq(1,length(Sorted_resids[[1]])),3)))))

rm(baseline_dummy)
names(Sorted_resis_points) = c("AUC_OD","drug","background","name")
Sorted_resis_points$AUC_OD = as.numeric(as.matrix(Sorted_resis_points$AUC_OD))

p4 <- ggplot(data = Sorted_resis_points, aes(x=drug, y=AUC_OD, fill=background)) +
  geom_bar(stat="summary",position=position_dodge(width=1)) +
  geom_point(aes(x=drug, y=AUC_OD, color=background), position = position_dodge(width = 1)) +
  scale_fill_manual(values=c(cols_light[1],'grey60')) +
  scale_color_manual(values=c(cols_main[1],'grey30'))

p4 + ggtitle(paste0("Screen hits, cutoff =",cutoff)) + ylim(0, 0.5)

#====== plot selected drugs OD curves ==========

sel_drugs = row.names(Sorted_resids)

drugs = substr(sel_drugs,1,(nchar(sel_drugs)-1))
concs = substr(sel_drugs,nchar(sel_drugs),nchar(sel_drugs))
wells = Map$Well[match(sel_drugs,paste0(Map$Drug,Map$ConcMock))]
Time_data = lapply(OD_data_corrected, `[[`, 1)
Time_data = as.data.frame(do.call(cbind, Time_data))

#Get OD for all wells across all plates
for(w in 1:length(wells))
{
  well = wells[w]
  
  well_i = match(well,names(OD_data_corrected[[1]]))
  well_data = lapply(OD_data_corrected, `[[`, well_i)
  well_data = as.data.frame(do.call(cbind, well_data))
  if(w==1)
  {Wells_data = list()}
  Wells_data[[length(Wells_data)+1]] = well_data
}
names(Wells_data) = wells

#Plot per well
par(mfrow=c(2,2),mar=c(4,4,4,4))
for(w in 1:length(wells))
{
  well = wells[w]
  drug = drugs[w]
  conc= concs[w]
  well_data = Wells_data[[w]]
  
  #Calculate error bars
  cbass_mean = as.data.frame(sapply(as.data.frame(t(well_data[1:3])),FUN=mean,na.rm=T))
  names(cbass_mean) = "CBASS_mean"
  cbass_sd = as.data.frame(sapply(as.data.frame(t(well_data[1:3])),FUN=sd,na.rm=T))
  names(cbass_sd) = "CBASS_sd"
  cbass_sd[1,1] = 0.001 #just to overcome a plotting error
  
  wt_mean = as.data.frame(sapply(as.data.frame(t(well_data[4:6])),FUN=mean,na.rm=T))
  names(wt_mean) = "wt_mean"
  wt_sd = as.data.frame(sapply(as.data.frame(t(well_data[4:6])),FUN=sd,na.rm=T))
  names(wt_sd) = "wt_sd"
  wt_sd[1,1] = 0.001 #just to overcome a plotting error
  
  highlight_col = cols_main[1]
  highlight_col_light = alpha(highlight_col,alpha=0.7)
  wt_col = "grey50"
  
  #Start plotting
  highlight_cols = c(rep(highlight_col,3),rep(wt_col,3))
  highlight_cols_light = c(rep(alpha(highlight_col,alpha=0.6),3),
                           rep(alpha(wt_col,alpha=0.6),3))
  
  here_pch = rep(c(15,19,17),2)
  here_pch_out = rep(c(0,1,2),2)
  
  #plot corrected OD values
  plot(Time_data[[1]],well_data[[1]],pch=here_pch[1],col=highlight_cols_light[1],
       ylab="OD595nm",xlab="Time (h)",main=paste0(drug," ",conc),type="p",ylim=c(0,0.8),cex=0.7)
  points(Time_data[[1]],well_data[[1]],pch=here_pch_out[1],col=highlight_cols[1],cex=0.7)
  for(i in 2:length(well_data))
  {
    points(Time_data[[i]],well_data[[i]],pch=here_pch[i],col=highlight_cols_light[i],type="p",cex=0.7)
    points(Time_data[[i]],well_data[[i]],pch=here_pch_out[i],col=highlight_cols[i],cex=0.7)
  }
  legend("topright",col=c(highlight_col,wt_col),legend = c("CBASS","WT"),pch=19,bty="n",cex=0.7)
  legend("topleft",col=wt_col,legend = c(1:3),pch=here_pch_out[1:3],bty="n",title = "Replicate",cex=0.7)
  
  #Plot error bars
  points(Time_data[[1]],cbass_mean[[1]],type="l",col=highlight_col,lwd=1)
  arrows(x0=Time_data[[1]],y0=(cbass_mean[[1]]-cbass_sd[[1]]),
         y1=(cbass_mean[[1]]+cbass_sd[[1]]),
         length=0.02, angle=90, code=3,lwd=0.5)
  
  points(Time_data[[4]],wt_mean[[1]],type="l",col=wt_col,lwd=1)
  arrows(x0=Time_data[[4]],y0=(wt_mean[[1]]-wt_sd[[1]]),
         y1=(wt_mean[[1]]+wt_sd[[1]]),
         length=0.02, angle=90, code=3,lwd=0.5)
}


dev.off()



