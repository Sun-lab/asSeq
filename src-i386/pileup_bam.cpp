/**********************************************************************
 *
 * asCounts_sam.c
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2011, Vasyl Zhabotynsky, UNC-CH
 *
 * first written Jul 11, 2011
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

void pileup_bam(char **Rinput,int* Rnum_inp,
				  char **Routput, int *Rmaxsize)
{
  /**
   *
   flag2kp  INT 	flags to keep
   flag2rm  INT 	flags to remove
   avgQ INT 	Minimum mapping quality allowed for a read to be used [0] 
   */

  int num_inp = *Rnum_inp;
  int maxsize = *Rmaxsize;
  int maxsize2 = maxsize*2;
  Rprintf("maxsize %d %d\n",maxsize,maxsize2);
  int currentpos,difmin;
  int currentmin=0;
  int flagmin=1,newmin=0;
  int* pileup_window;
  pileup_window = new int[maxsize2];
  for(int i=0;i<maxsize2;i++){
    pileup_window[i]=0;
  }
  int seqPosj, cigarPosj, lenj, offset1[MAX_LEN];
  char seqData[MAX_LEN];

  int isAccounted=0;
  
  int i, j,k,b,avg, flag, mapQ=0, keepIt1=1,chri=0;
  double favg;

  char *output;
  output  = Routput[0];
  string outputr,outputs,outputhap1,outputhap2,outputhap3;

  FILE *Fsnp, *FO1, *FO2; 
  //FILE *Fdump,*Fdump2;
  string dump,dump2;
  dump2=dump=string(output);
  dump+="_dump.bed";
  dump2+="_dump2.bed";
  //Fdump=fopen(dump.c_str(),"w");
  //Fdump2=fopen(dump2.c_str(),"w");
  //assert(Fdump);
  //assert(Fdump2);

  
  BamMultiReader reader;
  
  string currentchr="chr1";
  char str_buf[MAX_LEN], line[MAX_LEN],line2[MAX_LEN];  
  char chr[128], seq[255], qua[255];
  const char* strand="+-";
  //strand[0]='+';
  //strand[1]='-';

  std::vector<std::string> input,output_temp2,chrs;
  std:string temp;
  for(i=0;i<num_inp;i++){
    temp=Rinput[i];
    input.push_back(temp);
  }

  outputr=outputs=outputhap1=outputhap2=outputhap3=string(output);
  //outputr+="_reads.bed";
  outputs+="_pileup.bed";
  FO1=fopen(outputs.c_str(),"w");//to c string
  //FO2=fopen(outputr.c_str(),"w");//to c string
  assert(FO1);
  //assert(FO2);

  /*
   * M: alignment match (can be a sequence match or mismatch)
   *
   * I: insertion to the reference , e.g., 1I means 1 base pair
   *    in the sequence, but not in the reference genome
   *
   * D: deletion from the reference, e.g, 2D means 2 base pairs 
   *    in the refrence genome, but not in the sequence
   *
   * N: skipped region from the reference
  */  
  char *key, *alleles;
  key = (char*)calloc(255, 1);
  alleles = (char*)calloc(7, 1);

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
  /*
   * For a spit read, with starting psoition p0, 
   * the j-th position in the read may not be p0 + j - 1
   * and we specify it by p0 + offset1[j]
   *
   * well... here I assume the read length is smaller than 1023 bp
   */

  /**
   * read in the snp information
   */
  //Rprintf("snpList is %s\n", snpList);

  BamAlignment al1;
  vector<CigarOp>::iterator cigarIter;
  vector<CigarOp>::iterator cigarEnd;
  i=0;
  while(reader.GetNextAlignment(al1)){
	//if(i>1000000){break;}
    if((i++%1000000)==0)Rprintf("processing %d'th read\n",i);
    for (j=0; j<al1.Length; j++) { offset1[j] = j+1; }
    /* seqPosj is the current position in the sequence */
    seqPosj=0;    
    cigarIter = al1.CigarData.begin();
    cigarEnd  = al1.CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
      const CigarOp& op = (*cigarIter);
      switch (op.Type) {
      case (Constants::BAM_CIGAR_MATCH_CHAR)    :
      case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
      case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
        seqPosj += op.Length;
	    break;
      case (Constants::BAM_CIGAR_DEL_CHAR) :
      case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
        for (j=seqPosj; j<al1.Length; j++) { 
          offset1[j] += op.Length; 
        }
		break;
      case (Constants::BAM_CIGAR_INS_CHAR)      :
        for (j=seqPosj; j<(seqPosj+op.Length); j++) { 
          offset1[j] = -1;
        }
        for(j=(seqPosj+op.Length);j<al1.Length;j++){
          offset1[j]-= op.Length;
        }
        seqPosj  += op.Length;
        break;
        // shouldn't get here
      default:
		Rprintf("%c, %d", op.Type, op.Length);
        error("sorry, did not expect to see this character in cigar\n");
      }
    }
    newmin=al1.Position;
    difmin=newmin-currentmin;
	if((difmin<0)&(currentchr.compare(chrs[al1.RefID])==0)){
	  Rprintf("file should be sorted by position\n");
	  break;
	}
	//if((i%1000)==0){
	  //Rprintf("processing %d",i);
    //}
	//if(i==52095){Rprintf("%s %s %d %d\n",currentchr.c_str(),chrs[al1.RefID].c_str(),currentmin,difmin);}
    if((currentchr.compare(chrs[al1.RefID])!=0)|(difmin>maxsize)){	  	  
      if((currentchr.compare(chrs[al1.RefID])!=0)|(difmin>=maxsize2)){
	    //if(i==52095){Rprintf("got to block 1\n");}
        for(j=0;j<maxsize2;j++){
          if(pileup_window[j]>0){
            sprintf(line, "%s\t%d\t%d\n",currentchr.c_str(),currentmin+j,pileup_window[j]);
            fputs(line, FO1);
          }
          pileup_window[j]=0;
        }
      }else{
		//if(i==52095){Rprintf("got to block 2\n");}
        for(j=0;j<difmin;j++){
          if(pileup_window[j]>0){
            sprintf(line, "%s\t%d\t%d\n",currentchr.c_str(),currentmin+j,pileup_window[j]);
            fputs(line, FO1);
          }
          pileup_window[j]=0;
        }
		//if(i==52095){Rprintf("got to block 3\n");}
        for(j=difmin;j<maxsize2;j++){
          pileup_window[j-difmin]=pileup_window[j];
		  pileup_window[j]=0;
        }
      }
      currentmin=newmin;
      currentchr=chrs[al1.RefID];
    }
	//if(i==52095){Rprintf("%s %d\n",currentchr.c_str(),currentmin);}
	//if(i==52095){Rprintf("got to block 4\n");}
    for(j=0;j<al1.Length;j++){
      if(offset1[j]!=-1){
        currentpos=al1.Position+offset1[j]-currentmin;
		//if(i==52095){Rprintf("%d %d\n",currentpos,maxsize2);}
        if(currentpos<maxsize2){
          pileup_window[currentpos]++;
        }else{
          Rprintf("got out of window range\n");
        }
      }
    }	
	//if((i%1000)==0){
	  //Rprintf(" done with %d\n",i);
    //}
  }
  //Rprintf("checkpoint\n");
  for(j=0;j<maxsize2;j++){
    if(pileup_window[j]>0){
      sprintf(line, "%s\t%d\t%d\n",currentchr.c_str(),currentmin+j,pileup_window[j]);
      fputs(line, FO1);
    }
    pileup_window[j]=0;
  }
  
  reader.Close();
  Rprintf("%d lines read\n",i);
  fclose(FO1);
  //fclose(FO2);
  //fclose(Fdump);
  //fclose(Fdump2);

  delete [] pileup_window;  // When done, free memory pointed to by a.
  pileup_window = NULL; 
  free(key);
  free(alleles);
  Rprintf("total %d lines processed\n",i);
}
}// extern "C"
