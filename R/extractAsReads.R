`extractAsReads` <- 
function(input, snpList, outputTag=NULL, flag2kp=0, flag2rm=3844, 
         prop.cut=.5, min.avgQ=10, min.snpQ=10, phred=33, skip=0)
{
  if(any(is.null(input))|is.null(snpList)){
    stop("cannot continue without input,output and snpList files\n")
  }
  
  input = input[1]
  
  if(is.null(outputTag)){
    outputTag = strsplit(input[1],split="/")[[1]]
    outlen    = length(outputTag)
    outputTag[outlen] = strsplit(outputTag[outlen],split=".bam")[[1]]
    outputTag = paste(outputTag,collapse="/")
  }
  
  if(!all(file.exists(input))){ 
    stop("cannot find some of the input files: ", input[!file.exists(input)], "\n") 
  }
  
  if(!file.exists(snpList)){ 
    stop("cannot find file", snpList, "\n") 
  }
    
  Z = .C("extract_asReads", as.character(input),  as.character(outputTag), 
         as.character(snpList), as.double(prop.cut), as.double(min.avgQ), 
         as.double(min.snpQ),  as.integer(flag2kp), as.integer(flag2rm), 
         as.integer(phred), as.integer(skip), PACKAGE="asSeq")
  
  return(c(paste(outputTag,"hap1.bam",sep="_"),
           paste(outputTag,"hap2.bam",sep="_"),
           paste(outputTag,"hapN.bam",sep="_"),
           input))
}
