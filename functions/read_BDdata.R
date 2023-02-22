read_BDdata=function(path){
  directory=dir(path)
  list=list()
  for (i in 1:length(directory)) {
    dir_each=dir(paste0(path,"/",directory[i]))
    sample_name=paste0("s",str_split_fixed(directory[i],"_original",2)[1])
    matrix_name=dir_each[grep(pattern = "MolsPerCell",dir_each)]
    st_name=dir_each[grep(pattern = "Sample_Tag_Calls",dir_each)]
    data_matrix_i=read.csv(paste0(path,"/",directory[i],"/",matrix_name),
                           skip = 8,row.names = 1,header = T)
    sample_tag=read.csv(paste0(path,"/",directory[i],"/",st_name),
                        skip = 8,row.names = 1,header = T)
    data_matrix_i$sample_tag=sample_tag$Sample_Name
    data_matrix_i$sample_name=sample_name
    data_matrix_i_clean=subset(data_matrix_i,subset=!sample_tag%in%c("Multiplet","Undetermined"))
    row.names(data_matrix_i_clean)=paste0(sample_name,"-",row.names(data_matrix_i_clean))
    print(paste0("sample",i))
    multiplet_rate=round(as.numeric(table(sample_tag$Sample_Name)/nrow(sample_tag))[1],4)
    print(paste0("multiplet of ",sample_name," is ",multiplet_rate))
    print(paste0("effective cell rate is ",nrow(data_matrix_i_clean),
                 "/",nrow(data_matrix_i),"=",
                 round(nrow(data_matrix_i_clean)/nrow(data_matrix_i),4)))
    list=c(list(data_matrix_i_clean),list)
    names(list)[1]=sample_name
  }
  return(list)
}