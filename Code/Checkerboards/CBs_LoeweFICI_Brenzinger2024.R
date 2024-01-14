require(graphics)
require(stats)
require(Hmisc)
require(zoo)

#============ Specific_helper functions ===========
get_conc_grad <- function(c,steps=7)
{
  step = c/steps
  v=c((steps+1):1)
  v=v*step
  v = c-v
  v = round(c(v[2:(steps+1)],c),5)
  return(v)
}

get_CB_contours <- function(plot_data,Isoboles,drug_hrz,drug_vert,no_growth_limit = 0.000005)
{
  #remove non-growing concs
  if(length(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit))>0)
  {plot_data = plot_data[1:(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit)[1])]}
  plot_data = as.data.frame(t(plot_data))
  if(length(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit))>0)
  {plot_data = plot_data[1:(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit)[1])]}
  #plot_data = as.data.frame(t(plot_data))
  
  ctlns <- contourLines(x = as.numeric(row.names(plot_data)),
                        y = as.numeric(names(plot_data)),
                        z = as.matrix(plot_data), levels=Isoboles)
  return(ctlns)
}

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

get_feature <- function(plate_ids,feature,Plate_database)
{
  All_features = names(Plate_database)
  features = c()
  if(length(grep_exact(feature,All_features))>0)
  {
    column = grep_exact(feature,All_features)
    for(p in 1:length(plate_ids))
    {
      plate_id = plate_ids[p]
      features = c(features,as.character(Plate_database[grep_exact(plate_id,Plate_database[[1]]),column]))
    }
  }
  
  return(features)
}

#============ Set directories ===========
#Add local directory
here_path = ""

Load_dir = paste0(here_path,"")
Out_dir = paste0(here_path,"")

file_id = paste0(here_path,"PlateDatabase_combination.txt")
PlateDatabase = read.table (file_id, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

Plates = as.character(as.matrix(PlateDatabase$CB_id))
time_cutoffs = rep(10,length(Plates))
cb_wells = paste0(rep(LETTERS[1:8], each=12), rep(c(1:12),8))

#========== Get single drug MICs from file and normalize to MIC95 =======

file_id = paste0(here_path,"MIC_curves_combination.txt")
Drugs_MIC = read.table (file_id, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

strain_drugs = unique(paste0(Drugs_MIC$Strain,"_",Drugs_MIC$Drug))
for(i in 1:length(strain_drugs))
{
  strain_drug = strain_drugs[i]
  strain_drug_data = Drugs_MIC[grep_exact(strain_drug,paste0(Drugs_MIC$Strain,"_",Drugs_MIC$Drug)),]
  strain_drug_data = strain_drug_data[c(1,2)]
  names(strain_drug_data) = c("conc","fitness")
  MIC95 = round(strain_drug_data$conc[which.min(abs(strain_drug_data$fitness-0.05))],3)
  
  #normalize curve to MIC95
  strain_drug_data$conc = strain_drug_data$conc/MIC95
  
  if(i==1)
  {
    temp_Drugs_MIC = list()
    MIC95s = c()
  }
  temp_Drugs_MIC[[length(temp_Drugs_MIC)+1]] = strain_drug_data
  MIC95s = c(MIC95s,MIC95)
  
}
names(temp_Drugs_MIC) = strain_drugs
Drugs_MIC = temp_Drugs_MIC

MIC95s = as.data.frame(MIC95s)
row.names(MIC95s) = strain_drugs

rm(temp_Drugs_MIC,strain_drugs)

#========== Calculate grow (AUC)10h per plate: Growth ========

for(p in 1:length(Plates))
{
  plate_id = Plates[p]
  
  file_id = paste0(Load_dir,plate_id,".txt")
  dataset = read.table(file_id, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  time_cutoff = time_cutoffs[p]
  
  #Make 2nd point = average of 1 & 3 - noise because of the membrane bulging
  dataset[2,2:length(dataset)]=(dataset[1,2:length(dataset)]+dataset[3,2:length(dataset)])/2
  
  #Correct background by subtracting the 1st value
  corr_dataset=as.data.frame(t(as.data.frame(t(dataset[2:length(dataset)]))-t(dataset[3,2:length(dataset)])))
  corr_dataset=as.data.frame(cbind(dataset[[1]],corr_dataset))
  names(corr_dataset)[1] = "Time"
  dataset = corr_dataset
  rm(corr_dataset)
  
  growth = Trape_rule (as.data.frame(dataset[1:(length(grep(T,dataset[[1]]<time_cutoff))),]))
  growth = as.data.frame(growth)
  
  names(growth) = plate_id
  row.names(growth) = names(dataset)[2:length(dataset)]
  
  if(length(dataset)==97)
  {
    growth_96 = growth
    
    growth_384 = as.data.frame(rep(NA,384))
    row.names(growth_384)=paste0(rep(LETTERS[1:16],each=24),rep(seq(1,24),16),collapse = NULL)
    names(growth_384) = plate_id
    
    growth_384[match(row.names(growth_96),row.names(growth_384)),1] = growth_96
    
    growth = growth_384
    rm(growth_384,growth_96)
  }
  
  if(p==1)
  {Growth=growth} else
  {Growth = as.data.frame (cbind (Growth,growth))}
  
}
Growth[Growth<0]=0 #remove little noise brought in by slightly negative areas (consequence of background correction)

#============= Gather CBs data ===========

CBs = as.character(unique(PlateDatabase$CB_id))
bugs = get_feature(CBs,feature = "bug",Plate_database = PlateDatabase)
unique_bugs = unique(get_feature(CBs,feature = "bug",Plate_database = PlateDatabase))
CBs = CBs[order(match(bugs,unique_bugs))] 
rm(bugs,unique_bugs)

for(cb in 1:length(CBs))
{
  cb_id = CBs[cb]
  cb_rep = get_feature(cb_id,feature = "Replicate",PlateDatabase)
  
  bug=as.character(get_feature(cb_id,feature="bug",PlateDatabase))
  
  growth = Growth[match(cb_wells,row.names(Growth)),match(cb_id,names(Growth))]
  
  drug_hrz = get_feature(cb_id,feature="Horiz_drug",PlateDatabase)
  drug_vert = get_feature(cb_id,feature="Vert_drug",PlateDatabase)
  conc_drug_hrz = as.numeric(get_feature(cb_id,feature="Horiz_startC",PlateDatabase))
  conc_drug_vert = as.numeric(get_feature(cb_id,feature="Vert_startC",PlateDatabase))
  
  conc_drug_hrz = get_conc_grad(conc_drug_hrz,11)
  conc_drug_vert = get_conc_grad(conc_drug_vert,7)
  
  data = as.data.frame(t(matrix(growth,ncol=8,nrow=12)))
  data = data [rev(seq(1,12))]
  row.names(data) = rev(round(conc_drug_vert,3))
  names(data) = round(conc_drug_hrz,3)
  
  #calculate fitness
  data = data/data[8,1]
  
  #Store data
  if(cb==1)
  {CB_Growth = list()}
  
  CB_Growth[[length(CB_Growth)+1]] = data
}
names(CB_Growth) = CBs

#write Growth_vector
for(i in 1: length(CB_Growth))
{
  cb_id = names(CB_Growth)[i]
  cb_growth = as.data.frame(CB_Growth[[i]])
  cb_growth = rev(cb_growth)
  
  cb_growth_vector = as.numeric(as.matrix(t(cb_growth)))
  cb_conc_vector = paste0(rep(names(cb_growth),8),"_",
                         rep(row.names(cb_growth),each=12))
  cb_id_vector = rep(cb_id,length(cb_conc_vector))
  cb_wells = paste0(rep(LETTERS[1:8],each=12),rep(c(1:12),8))
  
  #get AUC from Growth, as CB_Growth is already normalized to no drug (fitness)
  growth = Growth[match(cb_wells,row.names(Growth)),match(cb_id,names(Growth))]
  
  cb_growth = as.data.frame(cbind(cb_id_vector,cb_wells,cb_conc_vector,growth))
  names(cb_growth) = c("cb_id","wells","Drugs conc (µg/ml), SMX_TMP","AUC10h")
  
  if(i==1)
    {Growth_vector = cb_growth} else
    {Growth_vector = as.data.frame(rbind(Growth_vector,cb_growth))}
}

file_id = paste0(Out_dir,"Growth_vector.txt")
write.table(Growth_vector,file=file_id,sep="\t",quote=F,row.names = F)
rm(file_id)

#============ Calculate & plot isobole contours =============

Isoboles= seq(0.1,0.5,by=0.05)
for(cb in 1:length(CBs))
{
  cb_id = CBs[cb]
  cb_rep = get_feature(cb_id,feature = "Replicate",PlateDatabase)
  
  bug=as.character(get_feature(cb_id,feature="bug",PlateDatabase))
  color = as.character(get_feature(cb_id,feature="Color",PlateDatabase))
  
  drug_hrz = get_feature(cb_id,feature="Horiz_drug",PlateDatabase)
  drug_vert = get_feature(cb_id,feature="Vert_drug",PlateDatabase)
  conc_drug_hrz = as.numeric(get_feature(cb_id,feature="Horiz_startC",PlateDatabase))
  conc_drug_vert = as.numeric(get_feature(cb_id,feature="Vert_startC",PlateDatabase))
  
  #get fitness
  data = CB_Growth[[grep_exact(cb_id,names(CB_Growth))]]
  
  #Get contours
  plot_data = data
  plot_data = plot_data[rev(1:length(plot_data[[1]])),]
  #CB_title = paste0(bug," rep",cb_rep," Bliss interaction")
  interac_ctlns = get_CB_contours(plot_data,Isoboles,drug_hrz,drug_vert)
  names(interac_ctlns) = Isoboles
  
  #normalize concs by MIC95 of wt
  MIC95s_st = MIC95s[grep(bug,row.names(MIC95s)),]
  MIC95s_st = as.data.frame(MIC95s_st)
  row.names(MIC95s_st) = row.names(MIC95s)[grep(bug,row.names(MIC95s))]
  MIC95_hrz = MIC95s_st[grep(drug_hrz,row.names(MIC95s_st)),1]
  MIC95_vert = MIC95s_st[grep(drug_vert,row.names(MIC95s_st)),1]
  
  for(i in 1:length(interac_ctlns))
  {
    interac_ctlns[[i]]$x = interac_ctlns[[i]]$x/MIC95_hrz
    interac_ctlns[[i]]$y = interac_ctlns[[i]]$y/MIC95_vert
  }
  
  if(cb==1)
  {All_interac_ctlns = list()}
  
  All_interac_ctlns[[length(All_interac_ctlns)+1]] = interac_ctlns
}
names(All_interac_ctlns) = CBs

#Write contours for source data
for(j in 1:length(All_interac_ctlns))
{
  interac_ctlns = All_interac_ctlns[[j]]
  for(i in 1:length(interac_ctlns))
  {
    isobole_data = interac_ctlns[[i]]
    isobole_SMX = isobole_data$x
    isobole_TMP = isobole_data$y
    isobole = rep(isobole_data$level,length(isobole_SMX))
    replicate = rep(names(All_interac_ctlns)[j],length(isobole_SMX))
    
    isobole_data = as.data.frame(cbind(replicate,isobole,isobole_SMX,isobole_TMP))
    
    if(j==1 && i==1)
    {Parsed_contours = isobole_data} else
    {Parsed_contours = as.data.frame(rbind(Parsed_contours,isobole_data))}
  }
}
names(Parsed_contours) = c("Strain_replicate","Isobole","SMX conc/SMX MIC","TMP conc/TMP MIC")

file_id = paste0(Out_dir,"Checkerboard countour lines.txt")
write.table(Parsed_contours,file=file_id,quote=F,row.names = F,sep="\t")
rm(file_id)

#========= Calculate FICI curves for isobole=0.1 (IC90) ===========
isobole=0.1
for(c in 1:length(All_interac_ctlns))
{
  cb = names(All_interac_ctlns)[c]
  
  data = All_interac_ctlns[[c]]
  bug = PlateDatabase$bug[match(cb,PlateDatabase$CB_id)]
  drug1 = PlateDatabase$Horiz_drug[match(cb,PlateDatabase$CB_id)]
  drug2 = PlateDatabase$Vert_drug[match(cb,PlateDatabase$CB_id)]
  replicate = PlateDatabase$Replicate[match(cb,PlateDatabase$CB_id)]
  
  data_isobole = data[[match(isobole,names(data))]]
  
  #Get MIC at isobole
  bugs_Drugs_MIC = unlist(strsplit(names(Drugs_MIC),split = "_"))
  bugs_Drugs_MIC = bugs_Drugs_MIC[seq(1,length(bugs_Drugs_MIC),by=2)]
  MIC_drug2 = names(Drugs_MIC)[grep_exact(bug,bugs_Drugs_MIC)]
  MIC_drug2 = MIC_drug2[grep(drug2,MIC_drug2)]
  MIC_drug1 = names(Drugs_MIC)[grep_exact(bug,bugs_Drugs_MIC)]
  MIC_drug1 = MIC_drug1[grep(drug1,MIC_drug1)]
  rm(bugs_Drugs_MIC)
  
  MIC_drug1 = Drugs_MIC[[match(MIC_drug1,names(Drugs_MIC))]]
  MIC_drug2 = Drugs_MIC[[match(MIC_drug2,names(Drugs_MIC))]]
  
  #get the closest fitness value to isobole
  drug1_alone = MIC_drug1$conc[which.min(abs(MIC_drug1$fitness - isobole))]
  drug2_alone = MIC_drug2$conc[which.min(abs(MIC_drug2$fitness - isobole))]
  
  #FICI_curve
  FICI_curve = data_isobole$x/drug1_alone + data_isobole$y/drug2_alone
  
  fici = min(FICI_curve)
  drug1_fici = data_isobole$x[which.min(FICI_curve)]
  drug2_fici = data_isobole$y[which.min(FICI_curve)]
  
  data_isobole = as.data.frame(cbind(FICI_curve,data_isobole$x,data_isobole$y,
                                     rep(cb,length(data_isobole$x)),
                                     rep(isobole,length(data_isobole$x))
  ))
  names(data_isobole) = c("FICi","SMX conc/SMX MIC","TMP conc/TMP MIC","Strain_replicate","isobole")
  
  if(c==1)
  {Data_isobole=data_isobole} else
  {Data_isobole=as.data.frame(rbind(Data_isobole,data_isobole))}
  
}

file_id = paste0(Out_dir,"FICi_data.txt")
write.table(Data_isobole,file = file_id,sep="\t",quote=F,row.names=F)
rm(file_id)



