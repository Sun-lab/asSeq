/**********************************************************************
 *
 * filter_bam.cc
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2011, Vasyl Zhabotynsky, UNC-CH
 *
 * first written Jul 11, 2011
 * 
 * modified by Chong Jin 2016-06-28
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 **********************************************************************/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include "stdhash.hh"
#include "BamMultiReader.h"
#include "BamWriter.h"
#include "BamConstants.h"
#include "bamtools_sort.h"
using namespace BamTools;
#include <iostream>
#include <vector>
#include <string>
using namespace std;
#include "alName.h"

#define MAX_LEN 2048

extern "C" {
/**********************************************************************
 *
 * generate counts from the ouptut of tophat: accepted_hits.sam
 *
 **********************************************************************/

/**
 *
 
 Three lines of the input file looks like the following:
 
 IL26_1382:1:109:660:632 137     chr1    309     1       36M     *       0       0
 CAACCCCAACCCTAACCCTAACCCCTAACCCTAACC    <<=<<<<<<<4<<<9<;<<9<<;<;<;7<;<9:0:<
 NM:i:1  NH:i:4  CC:Z:chr2       CP:i:114076909
 
 IL26_1382:1:61:717:749  99      chr1    3037    3       36M     =       3093    0
 TGGGGAGGCAGCTGTAACTCAAAGCCTTAGCCTCTG    >>>>><>>>3>>;>:96=>9>2>><:8<39:+1::5
 NM:i:0  NH:i:2  CC:Z:chr9       CP:i:3287
 
 IL26_1382:1:61:717:749  147     chr1    3093    3       36M     =       3037    0
 TCAGGCGCCAAAGGGATTCTGCCAGCATAGTGCTCC    3;;&<6.><<<>><>7<<<<><>>>><<<>><<>>>
 NM:i:1  NH:i:2  CC:Z:chr9       CP:i:3343
 
 Each line of maq mapview output consists of
 [1] QNAME
 [2] FLAG
 [3] RNAME   reference name (chromsome)
 [4] POS     leftmost position of the sequence
 [5] MAPQ    phred-scaled quality score
 [6] CIGAR   M: match
 [7] MRNM    mate reference sequence
 [8] MPOS    left most position of mate sequence
 [9] ISIZE   inferred insert size
 [10] SEQ    Sequence
 [11] QUAL   ASCII - 33 gives Phred base quality
 [12] TAG:VTYPE:VALUE
             (1) NM: number of nucleotide difference
             (2) NH: number of alignments
             (3) CC: reference name of next hit
             (4) CP: Leftmost coordinate of the next hit
 */

/* function is defined in asCounts2_bam.cc */
int StringToNumber ( const string &Text );

void filter_bam(char **Rinput,int *Rnum_inp, char **Routput, 
                double *Rmin_avgQ, int *Rmin_mapQ, 
                int *Rflag2kp, int *Rflag2rm, int *Rphred)
{
  /**
   *
   flag2kp  INT 	flags to keep
   flag2rm  INT 	flags to remove
   min_avgQ INT 	Minimum mapping quality allowed for a read to be used [0] 
   */
  int delim_ind=0;//this addition
  int phred=*Rphred;
  Rprintf("phred correction is %d\n",phred);
  
  int min_mapQ = *Rmin_mapQ;
  int flag2kp  = *Rflag2kp;
  int flag2rm  = *Rflag2rm;
  int num_inp  = *Rnum_inp;
  double min_avgQ = *Rmin_avgQ;
  double prop1,prop2;

  uint32_t tagval;
    
  int i, j,avg, mapQ=0, keepIt1=1,keepIt2=1,flag=0,chri=0;
  char line[MAX_LEN];
  double favg;
  int skip0=0; // the number of reads skipped due to flags
  int skip1=0; // the number of reads skipped due to min_avgQ
  int skip2=0; // the number of reads skipped due to min_mapQ
  int skip3=0; // skipped due to malformat
  int isSame=0;

  const char* strand="+-";
    
  char* pBaseQualities;
  char* output;
  
  //FILE *Finfo;
  output  = Routput[0];
  string outputf;
  outputf = string(output);

  //outputi+="_info.bed";
  //Finfo=fopen(outputi.c_str(),"w");
  //assert(Finfo);

  BamMultiReader reader;
  BamWriter writer;
  
  std::vector<std::string> input,chrs;
  std:string temp;
  for(i=0;i<num_inp;i++){
    temp=Rinput[i];
    input.push_back(temp);
  }

  if(!reader.Open(input)){
	  Rprintf("Could not open some input files.\n");
	  return;
  }
  
  const SamHeader header=reader.GetHeader();
  const RefVector references=reader.GetReferenceData();
  SamSequenceConstIterator begin = header.Sequences.ConstBegin();
  SamSequenceConstIterator iter  = begin;
  SamSequenceConstIterator end   = header.Sequences.ConstEnd();
  
  for ( ; iter != end; ++iter ) {
    const SamSequence& currentSeq = (*iter);
    chrs.push_back(currentSeq.Name);
    chri++;
  }
  
  if(!writer.Open(outputf,header,references)){
	  Rprintf("Could not open %s output file.\n",outputf.c_str());
	  return;
  }
  
  BamAlignment al1,al2;
  vector<CigarOp>::iterator cigarIter;
  vector<CigarOp>::iterator cigarEnd;
  i=0;
  
  while(reader.GetNextAlignment(al1)){
    i++;
    keepIt1=1;
    //cout << i<< endl;
    if(al1.Name.empty()){
      al1.SetIsPaired(false);
      skip3++;
      keepIt1=0;
      //cout << "skipping " << i << endl;
    }
    delim_ind=whichDelim(al1.Name);

    //if(i>100)break;
    if((i%1000000)==0)Rprintf(".");
    if((i%10000000)==0)Rprintf("\nprocessing %d million-th line\n", i/1000000);
    
    // keepIt1 and keepIt2 indicate whether to keep 1st or 2nd read
    
    if(al1.IsPaired()){	
      //cout << "is paired";  
      if(reader.GetNextAlignment(al2)){
        i++;
        //cout << i<< endl;
        if((i%1000000)==0)Rprintf(".");
        if((i%10000000)==0)Rprintf("\nprocessing %d million-th line\n", i/1000000);

        isSame=((split(al1.Name,DELIM_LOOKUP[delim_ind])[0]).compare(split(al2.Name,DELIM_LOOKUP[delim_ind])[0]) == 0);
        
        if (! isSame) {
          cout<<al1.Name<<endl;
          cout<<al2.Name<<endl;
          error("expect a paired-read, but see something else\n");
        }
        
        keepIt2=1;
      }else{
        keepIt2=0;
      }
    }else{
      keepIt2=0;
    }
    
    if ((keepIt1==1)&&(((al1.AlignmentFlag & flag2kp) != flag2kp) || (al1.AlignmentFlag & flag2rm))){
      skip0  += 1;
      keepIt1 = 0;
    }
    
    if ((keepIt2==1)&&(((al2.AlignmentFlag & flag2kp) != flag2kp) || (al2.AlignmentFlag & flag2rm))){
      skip0  += 1;
      keepIt2 = 0;
    }
    
    if(keepIt1==1){
      avg=0;
      pBaseQualities = (char*)al1.Qualities.data();
      
      for(j=0;j<al1.Length;j++){
        avg += pBaseQualities[j];
      }
      
      favg=(double)avg/(al1.Qualities).size()-phred;

      if(favg < min_avgQ){
        skip1  += 1;
        keepIt1 = 0;
      }
    }
    
    if(al1.MapQuality < min_mapQ){
      skip2  += 1;
      keepIt1 = 0;
    }

    if(keepIt2==1){
      avg=0;
      pBaseQualities = (char*)al2.Qualities.data();
      
      for(j=0;j<al2.Length;j++){
        avg += pBaseQualities[j];
      }
      
      favg=(double)avg/(al2.Qualities).size()-phred;

      if(favg < min_avgQ){
        skip1  += 1;
        keepIt2 = 0;
      }
    }
    
    if(al2.MapQuality < min_mapQ){
      skip2 += 1;
      keepIt2=0;
    }
    
    if((keepIt1==0) && (keepIt2==0)){
      continue;
    }else if((keepIt1!=0) && (keepIt2==0)){
      // "If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80." -- Sam format specification
      // Unset 0x8 to avoid "Mate unmapped flag should not be set for unpaired reads" error in picard.
      //al1.SetIsMateMapped(true);//~BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
      al1.InsertSize   = 0;
      al1.MateRefID    = -1;
      al1.MatePosition = -1;
      al1.SetIsPaired(false);//BAM_ALIGNMENT_PAIRED=0x0001	      
      writer.SaveAlignment(al1);
    }else if((keepIt1==0) && (keepIt2!=0)){
      //al2.SetIsMateMapped(true);//~BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
      al2.InsertSize   = 0;
      al2.MateRefID    = -1;
      al2.MatePosition = -1;
      al2.SetIsPaired(false);//BAM_ALIGNMENT_PAIRED=0x0001	      
      writer.SaveAlignment(al2);
    }else{
      //al1.SetIsMateMapped(true);//~BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
      //al2.SetIsMateMapped(true);//~BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
      writer.SaveAlignment(al1);
      writer.SaveAlignment(al2);
    }
  }
  
  reader.Close();
  writer.Close();
  
  Rprintf("%d reads are skipped due to flags\n",skip0);
  Rprintf("%d reads are skipped due to min_avgQ\n",skip1);
  Rprintf("%d reads are skipped due to min_mapQ\n",skip2);
  Rprintf("%d reads were malformatted\n",skip3);
  Rprintf("%d reads are processed in total\n",i);
  
}
}// extern "C"
