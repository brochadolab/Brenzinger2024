require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(ggplot2)
require(caTools)
require(pracma)
library("drc")

#=========== Helper functions =======
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

Trape_rule <- function(dataset)
{
  areas = c()
  
  xData = as.numeric(as.vector(dataset[[1]]))
  nr_points = length(xData)
  h = (xData[length(xData)]-xData[1])/(nr_points-1)
  
  for(i in 2:length(dataset))
  {
    yData = as.numeric(as.vector(dataset[[i]]))
    
    t = (yData[1] + yData[length(yData)])/2
    for(j in 2:(length(yData)-1))
    {
      t = t + yData[j]
    }
    
    Tn = t*h
    areas = c(areas,Tn)
  }
  names(areas) = names(dataset)[2:length(dataset)]
  return(areas)
  
}

#============ Set the stage ===========
here_path = "/"

file_id = paste0(here_path,"MasterMIC_database.txt")
MasterMIC_database = read.table (file_id, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
MasterMIC_database = MasterMIC_database[MasterMIC_database$Keep=="keep",]

Load_dir = paste0(here_path,"")
Out_dir = paste0(here_path,"")

file_id = paste0(Load_dir,"Parsed_MIC_data.txt")
MIC_data = read.table (file_id, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
#============ Calculate AUC, plot MIC curves and calculate MIC_target  ===========

#Get AUCs per MIC_id
All_mics = MasterMIC_database$MIC_id
targetMIC=0.90
MICs = c()
fit_MICs = c()

pdf(paste0(Out_dir,"MIC_curves_95.pdf"),useDingbats = F)
par(mfrow=c(4,4),mar=c(2,2,2,2))
for(m in 1:length(All_mics))
{
  mic_id = All_mics[m]
  mic_data = MIC_data[grep_exact(mic_id,MIC_data$MIC_id),]
  mic_time = mic_data$Time_h
  mic_OD = mic_data[grep("OD",names(mic_data))]
  mic_conc = mic_data[-c(grep("OD",names(mic_data)),grep("Time",names(mic_data)),
                         grep("MIC",names(mic_data)))]
  
  concs = round(mic_conc[1,],3)
  cols = colorRampPalette(c("#FFA500","white"))(length(concs)+2)
  plot(mic_time,mic_OD[[1]],type="b",pch=19,col=cols[1],
       ylab="OD",xlab="Time_h",ylim=c(0,1),main=mic_id)
  
  for(i in 2:length(mic_OD))
  {points(mic_time,mic_OD[[i]],type="b",pch=19,col=cols[i])}
  legend("topleft",legend = concs,pch=19,col=cols,bty="n")
  
  #mic_curve2 = as.data.frame(cbind(t(concs),sapply(mic_OD,FUN=trapz,x=mic_time)))
  mic_curve = as.data.frame(cbind(t(concs),
                                  Trape_rule(as.data.frame(cbind(mic_time,mic_OD)))))
  mic_curve$V2[mic_curve$V2<0]=0
  
  names(mic_curve) = c("conc","auc")
  mic_curve = as.data.frame(cbind(mic_curve,mic_curve$auc/max(mic_curve$auc,na.rm = T)))
  names(mic_curve) = c("conc","auc","fitness")
  
  mL <- drm(mic_curve$fitness ~ mic_curve$conc, fct = L.3(), type = "continuous")
  new_concs = seq(0,max(mic_curve$conc,na.rm = T)+200,by=0.001)
  predict_data = predict(mL,newdata = as.data.frame(new_concs,by=0.1))
  
  mic = new_concs[which.min(abs(predict_data-(1-targetMIC)))]
  
  upper_xlim = 1.2*(max(max(mic_curve$conc,na.rm = T),mic))
  #plot(mic_curve$conc,mic_curve$auc,pch=19,col=cols[1])
  plot(mic_curve$conc,mic_curve$fitness,pch=19,col=cols[1],
       ylim=c(0,1.0),xlim=c(0, upper_xlim), xlab="drug conc (mg/l)",ylab="fitness")
  points(new_concs,predict_data,type="l")
  points(new_concs[which.min(abs(predict_data-(1-targetMIC)))],
         predict_data[which.min(abs(predict_data-(1-targetMIC)))],pch=19)
  
  MICs = c(MICs,mic)
  fit_MICs = c(fit_MICs,predict_data[which.min(abs(predict_data-(1-targetMIC)))])
  
}
dev.off()

MasterMIC_database = as.data.frame(cbind(MasterMIC_database,MICs,fit_MICs))
names(MasterMIC_database)[c(length(MasterMIC_database)-1,length(MasterMIC_database))] = c("MIC","Fitness at MIC")
file_id = paste0(Out_dir,"Estimated_MIC_data.txt")
write.table(MasterMIC_database,file = file_id,quote=F,sep = "\t",row.names = F,col.names = T)
