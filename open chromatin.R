# find connection between promoter/enhancer from different samples and Epitensor interaction

# import and filter ATAC-seq data
# ATAC-seq data address: /share/resource3/data/TCGA/ACCx/ATACseq/4.PeakCalls_macs2_blacklistfiltered/
atac = read.table("/Users/peter/Desktop/UCSD_BISB/Fall_Rotation/4.PeakCalls_macs2_blacklistfiltered/fbd18a1c-d602-488f-a072-91c61919a171_BLfiltered.bed",
                   header=FALSE, sep="\t")
write.table(atac[,c(1:3,5)], file="/Users/peter/Desktop/UCSD_BISB/Fall_Rotation/promoter_enhancer/atac_fbd.bed", col.names=F,row.names=F,quote=F,sep='\t')

# Command to generate >=50% overlaps between promoter/enhancer region and Epitensor interaction regions
# input atac_fbd.bed: ATAC-seq data
# input hg38.loci1/2.bed: first and second loci from hg38 Epitensor pair data
# output overlap_fbd.bed: overlaped region between ATAC-seq data and first and second loci of Epitensor data

# bedtools intersect -wa -wb -a atac_fbd.bed -b hg38.loci1.bed hg38.loci2.bed -names 1st_loci 2nd_loci -f 0.5 > overlap_fbd.bed

# find ATAC-seq data paired in Epitensor
library(dplyr)
overlap = read.table("/Users/peter/Desktop/UCSD_BISB/Fall_Rotation/promoter_enhancer/overlap_fbd.bed",header=FALSE, sep="\t")
first_intersect = overlap[overlap[,5]=="1st_loci",9] # ATAC's overlap with first loci in Epitensor pairs
second_intersect = overlap[overlap[,5]=="2nd_loci",9] # ATAC's overlap with second loci in Epitensor pairs
inter = intersect(first_intersect, second_intersect)
atac_1 = overlap[overlap[,9] %in% inter,]
atac_1 = atac_1[order(atac_1$V9, decreasing=FALSE),]
atac_1 = atac_1[,c(1:5,9)]
atac_1[1:5,]

# find pairs of promoter and enhancer peaks and their corresponding TSS
tss = read.table("/Users/peter/Desktop/UCSD_BISB/Fall_Rotation/promoter_enhancer/valid_TSS.txt",header=T,sep='\t')
id = unique(atac_1[,6])
valid = rep(FALSE, length(id))
select_tss = c()
start_time = Sys.time()
for(i in 1:length(id)){
  fir = atac_1[atac_1[,6]==id[i] & atac_1[,5]=="1st_loci",]
  sec = atac_1[atac_1[,6]==id[i] & atac_1[,5]=="2nd_loci",]
  temp1 = tss[tss[,1]==fir[1,1] & tss[,5]=="+" & (tss[,2]-2000<=max(fir[,3]) | tss[,3]+2000>=min(fir[,2])),]
  temp2 = tss[tss[,1]==fir[1,1] & tss[,5]=="-" & (tss[,2]-2000<=max(sec[,3]) | tss[,3]+2000>=min(sec[,2])),]
  if(nrow(temp1)>0 & nrow(temp2)>0){
    pro_fir = rep(FALSE, nrow(temp1)+nrow(temp2))
    pro_sec = rep(FALSE, length(pro_fir))
    for(j in 1:nrow(fir)){
      a = c(!((temp1[,2]-fir[j,3])>2000|(fir[j,2]-temp1[,3])>500), !((temp2[,2]-fir[j,3])>500|(fir[j,2]-temp2[,3])>2000))
      pro_fir = pro_fir | a
    }
    for(j in 1:nrow(sec)){
      b = c(!((temp1[,2]-sec[j,3])>2000|(sec[j,2]-temp1[,3])>500), !((temp2[,2]-sec[j,3])>500|(sec[j,2]-temp2[,3])>2000))
      pro_sec = pro_sec | b
    }
    if(1 %in% c(pro_fir+pro_sec)){
      valid[i] = TRUE
      select_tss[[i]] = which(c(pro_fir+pro_sec) %in% 1)
    }
  }
  rec_time = Sys.time()
  print(paste("progress: ", i/length(id), "; remain time: ", (rec_time-start_time)*(length(id)-i)/i, sep=""))
}
valid
select_tss
valid_region = id[valid] # valid id of epitensor pairs
atac_2 = atac_1[atac_1[,6] %in% id[valid],] #combination of promoter and enhancer peaks, atac_2 only has only valid promoters and enhancers


# isolate promoter and enhancer peaks
library(dplyr)
atac_3 = atac_2[order(atac_2[,6], atac_2[,5]),]
is_pro = c()
is_enh = c()
start_time = Sys.time()
for(i in 1:length(valid_region)){ # length(valid_region)
  fir = atac_3[atac_3[,6]==valid_region[i] & atac_3[,5]=="1st_loci",]
  sec = atac_3[atac_3[,6]==valid_region[i] & atac_3[,5]=="2nd_loci",]
  temp1 = tss[tss[,1]==fir[1,1] & tss[,5]=="+",]
  temp2 = tss[tss[,1]==fir[1,1] & tss[,5]=="-",]
  temp_tss = rbind(temp1, temp2)
  temp_tss = temp_tss[select_tss[[which(id %in% valid_region[i])]],]
  temp1 = temp_tss[temp_tss[,5]=="+",]
  temp2 = temp_tss[temp_tss[,5]=="-",]
  
  pro_fir = rep(FALSE, nrow(temp_tss))
  pro_sec = rep(FALSE, nrow(temp_tss))
  for(j in 1:nrow(fir)){
    a = c(!((temp1[,2]-fir[j,3])>2000|(fir[j,2]-temp1[,3])>500), !((temp2[,2]-fir[j,3])>500|(fir[j,2]-temp2[,3])>2000))
    pro_fir = pro_fir | a
  }
  for(j in 1:nrow(sec)){
    b = c(!((temp1[,2]-sec[j,3])>2000|(sec[j,2]-temp1[,3])>500), !((temp2[,2]-sec[j,3])>500|(sec[j,2]-temp2[,3])>2000))
    pro_sec = pro_sec | b
  }
  
  fir_pro = rep(FALSE, nrow(fir))
  fir_enh = rep(FALSE, nrow(fir))
  sec_pro = rep(FALSE, nrow(sec))
  sec_enh = rep(FALSE, nrow(sec))
  
  for(j in 1:length(pro_fir)){
    if(pro_fir[j]){ # 1st loci: 2nd loci = promoter : enhancer; check if 1st loci are in promoter range of temp_tss[j,]
      if(temp_tss[j,5]=="+"){
        fir_pro = fir_pro | (!(temp_tss[j,2]-fir[,3]>2000|fir[,2]-temp_tss[j,3]>500))
        sec_enh = sec_enh | (temp_tss[j,2]-sec[,3]>2000|sec[,2]-temp_tss[j,3]>500)
      }
      else{
        fir_pro = fir_pro | (!(temp_tss[j,2]-fir[,3]>500|fir[,2]-temp_tss[j,3]>2000))
        sec_enh = sec_enh | (temp_tss[j,2]-sec[,3]>500|sec[,2]-temp_tss[j,3]>2000)
      }
    }
    else{ # 1st loci : 2nd loci = enhancer : promoter
      if(temp_tss[j,5]=="+"){
        sec_pro = sec_pro | (!(temp_tss[j,2]-sec[,3]>2000|sec[,2]-temp_tss[j,3]>500))
        fir_enh = fir_enh | (temp_tss[j,2]-fir[,3]>2000|fir[,2]-temp_tss[j,3]>500)
      }
      else{
        sec_pro = sec_pro | (!(temp_tss[j,2]-sec[,3]>500|sec[,2]-temp_tss[j,3]>2000))
        fir_enh = fir_enh | (temp_tss[j,2]-fir[,3]>500|fir[,2]-temp_tss[j,3]>2000)
      }
    }
  }
  is_pro = c(is_pro, c(fir_pro, sec_pro))
  is_enh = c(is_enh, c(fir_enh, sec_enh))
  rec_time = Sys.time()
  print(paste("progress: ", i/length(valid_region), "; remain time: ", (rec_time-start_time)*(length(valid_region)-i)/i, sep=""))
}
atac_3 = cbind(atac_3, is_pro, is_enh)
write.table(atac_3, file="/Users/peter/Desktop/UCSD_BISB/Fall_Rotation/promoter_enhancer/fbd_proenh_peaks_for_motif.bed",col.names=T,row.names=F,quote=F,sep='\t')

temp1 = tab3[atac_3[,7],1:6]
temp2 = tab3[atac_3[,8],1:6]
write.table(unique(temp1[,1:4]), file="/Users/peter/Desktop/UCSD_BISB/Fall_Rotation/promoter_enhancer/fbd_unique_promoter_peaks_for_motif.bed",col.names=F,row.names=F,quote=F,sep='\t')
write.table(unique(temp2[,1:4]), file="/Users/peter/Desktop/UCSD_BISB/Fall_Rotation/promoter_enhancer/fbd_unique_enhancer_peaks_for_motif.bed",col.names=F,row.names=F,quote=F,sep='\t')

