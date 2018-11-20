`prepareBAM` <-
function(input, outputTag, sortIt=TRUE, getUniqMapping=TRUE, fixmate=FALSE,
filterIt=TRUE, sort.nmax=NULL, sort.mem=NULL, min.avgQ=10, 
min.mapQ=10, flag2kp=0, flag2rm=1796, phred=33, 
keepSortedBAM=FALSE, keepUniqBAM=FALSE)
{
  
  if(length(input) > 1){
    warning("there are more than one input files, only the first one is used\n")
  }
  input = input[1]
  
  if(length(outputTag) > 1){
    warning("there are more than one outputTags, only the first one is used\n")
  }
  outputTag = outputTag[1]
  
  if(!is.null(sort.nmax) && !is.numeric(sort.nmax)){
    stop("sort.nmax should be either numeric or NULL\n")
  }
  
  if(!is.null(sort.mem) && !is.numeric(sort.mem)){
    stop("sort.mem should be either numeric or NULL\n")
  }
  
  # ------------------------------------------------------------------------
  # sort the bam file by name first 
  # ------------------------------------------------------------------------
  
  if(sortIt){
    outputSorted = sprintf("%s_sorted_by_name.bam", outputTag)
    
    params = c(getwd(),"sort","-in",input,"-out",outputSorted,"-byname")
    
    if(!is.null(sort.nmax)){
      params = c(params,"-n",as.character(sort.nmax))
    }
    
    if(!is.null(sort.mem)){
      params = c(params,"-mem",as.character(sort.mem))
    }
    
    nparams = length(params)
    
    Z = .C("sort_bam", as.character(params), as.integer(nparams), PACKAGE="asSeq")
  }
  
  # ------------------------------------------------------------------------
  # get uniquely mapped reads. Otherwise those reads mapped to multiple 
  # places in the genome will cause problems in the counting step
  # ------------------------------------------------------------------------
  
  if(getUniqMapping){
    
    if(sortIt){
      inputUniq  = outputSorted
      outputUniq = sprintf("%s_sorted_by_name_uniq.bam", outputTag)
    }else{
      inputUniq  = input
      outputUniq = sprintf("%s_uniq.bam", outputTag)
    }
    
    num_files = 1
    Z = .C("getUnique_bam", as.character(inputUniq), as.integer(num_files),
           as.character(outputUniq), as.integer(fixmate), PACKAGE="asSeq")
  }
  
  # ------------------------------------------------------------------------
  # get uniquely mapped reads. Otherwise those reads mapped to multiple 
  # places in the genome will cause problems in the counting step
  # ------------------------------------------------------------------------
  
  if(filterIt){
    
    if(sortIt && getUniqMapping){
      inputFilter  = outputUniq
      outputFilter = sprintf("%s_sorted_by_name_uniq_filtered.bam", outputTag)
    }else if((!sortIt) && getUniqMapping){
      inputFilter  = outputUniq
      outputFilter = sprintf("%s_uniq_filtered.bam", outputTag)
    }else if(sortIt && (!getUniqMapping)){
      inputFilter  = outputSorted
      outputFilter = sprintf("%s_sorted_by_name_filtered.bam", outputTag)
    }else{
      inputFilter  = input
      outputFilter = sprintf("%s_filtered.bam", outputTag)
    }
        
    num_files = 1
    
    Z = .C("filter_bam", as.character(inputFilter),as.integer(num_files),
           as.character(outputFilter), as.double(min.avgQ), 
           as.integer(min.mapQ),as.integer(flag2kp), as.integer(flag2rm), 
           as.integer(phred), PACKAGE="asSeq")
  }
    
  # ------------------------------------------------------------------------
  # get uniquely mapped reads. those reads mapped to multiple places in the 
  # genome is harder to handel in the counting step
  # ------------------------------------------------------------------------
    
  if(sortIt && (!keepSortedBAM)){
    system(sprintf("rm %s", outputSorted))
  }
  
  if(getUniqMapping && (!keepUniqBAM)){
    system(sprintf("rm %s", outputUniq))
  }
  
  return (1)
}
