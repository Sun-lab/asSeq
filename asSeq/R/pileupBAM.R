`pileupBAM` <-
function(input, output, maxsize=1E6){
  
  if(length(input) > 1){
    warning("there are more than one input files, only the first one is used\n")
  }
    
  input = input[1]
  
  if(!file.exists(input)){ 
    stop("cannot find the input file: ", input, "\n") 
  } 
  
  num_files = 1
  
  Z = .C("pileup_bam", as.character(input), as.integer(num_files), 
         as.character(output), as.integer(maxsize), PACKAGE="asSeq")
  
  return(c(input, output))
}
