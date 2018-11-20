/**********************************************************************
 *
 * extract_asReads.cc
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2011, Vasyl Zhabotynsky, UNC-CH
 *
 * first written Jul 11, 2011
 * modified by Wei Sun and Vasyl Zhabotynsky 2012-11-06
 * modified by Wei Sun 2014-02-13
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

hash_map_char<char*> *load_snp_set(FILE *fp)
{
  char chr[255], allele1[7], allele2[7];
  char buffer[2], key[255];
  int pos, c;
  hash_map_char<char*> *hash = new hash_map_char<char*>;
  while (fscanf(fp, "%s\t%d\t%s\t%s", chr, &pos, allele1, allele2) == 4) {
    
    /* skip it if this is not a heterzygous SNP*/
    if(strcmp(allele1, allele2)==0){ continue; }
    
    sprintf(key, "%s.%d", chr, pos);
    buffer[0] = allele1[0];
    buffer[1] = allele2[0];
    buffer[2] = 0;
    while ((c = fgetc(fp)) != EOF && c != '\n');
    hash->insert(key, strdup(buffer));
  }
  return hash;
}
  
void extract_asReads(char **Rinput, char **Routput, char ** RsnpList, 
                  double* Rpropcut,double* RavgQ, double* RsnpQ, 
                  int *Rflag2kp, int *Rflag2rm, int* Rphred, int *skip)
{
  /**
   *
   flag2kp  INT 	flags to keep
   flag2rm  INT 	flags to remove
   avgQ INT 	Minimum mapping quality allowed for a read to be used [0] 
   */
  
  // proportion of AS reads to haplotype 1 and 2.
  int delim_ind=0;//this addition
  double prop1, prop2;
  
  double avgQ    = *RavgQ;
  double snpQ    = *RsnpQ;
  double propcut = *Rpropcut;

  int phred      = *Rphred;
  int flag2kp    = *Rflag2kp;
  int flag2rm    = *Rflag2rm;
  
  long int pos1[MAX_LEN], pos2[MAX_LEN], pos3[MAX_LEN], currentpos;
  int seqPosj, cigarPosj, lenj, offset1[MAX_LEN], offset2[MAX_LEN];
  char seqData[MAX_LEN];

  int isAccounted=0;
  int isSame=0;
  
  int i, j,k,b,avg, flag, mapQ=0, keepIt1=1,keepIt2=1,k1,k2,k3,kall,chri=0;
  double favg;
	int nSNP =0; // total number of SNP harboring reads
	int skip0=0; // the number of reads skipped due to flags
	int skip1=0; // the number of reads skipped due to mapQ
	int skip2=0; // the number of reads skipped due to quality at SNP position
    
  char *output, *snpList;
  output  = Routput[0];
  string outputhap1, outputhap2, outputhap3;

  FILE *Fsnp; 
  
  BamMultiReader reader;
  BamWriter writer1,writer2,writer3;
  char* pBaseQualities;
  
  char currentchr[MAX_LEN];
  char str_buf[MAX_LEN], line[MAX_LEN],line2[MAX_LEN];  
  char chr[128], seq[255], qua[255];
  const char* strand="+-";
  //strand[0]='+';
  //strand[1]='-';

  std::vector<std::string> input,chrs;
  input.push_back(Rinput[0]);

  outputhap1=outputhap2=outputhap3=string(output);
  outputhap1+="_hap1.bam";
  outputhap2+="_hap2.bam";
  outputhap3+="_hapN.bam";
  
  snpList = RsnpList[0];

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
  key     = (char*)calloc(255, 1);
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
  
  if(!writer1.Open(outputhap1,header,references)){
	  Rprintf("Could not open %s output file.\n",outputhap1.c_str());
	  return;
  }
  
  if(!writer2.Open(outputhap2,header,references)){
	  Rprintf("Could not open %s output file.\n",outputhap2.c_str());
	  return;
  }
  
  if(!writer3.Open(outputhap3,header,references)){
	  Rprintf("Could not open %s output file.\n",outputhap3.c_str());
	  return;
  }
  
  int n_unrecognized = 0;
  
  /*
   * For a split read, with starting psoition p0, 
   * the j-th position in the read may not be p0 + j - 1
   * and we specify it by p0 + offset1[j]
   *
   * well... here I assume the read length is smaller than 1023 bp
   */

  /**
   * read in the snp information
   */
  
  //Rprintf("snpList is %s\n", snpList);

  Fsnp = fopen(snpList, "r");
  assert(Fsnp);
	hash_map_char<char*> *hash_map = 0;
  hash_map = load_snp_set(Fsnp);
  fclose(Fsnp);

  BamAlignment al1,al2;
  vector<CigarOp>::iterator cigarIter;
  vector<CigarOp>::iterator cigarEnd;
  
  i=0;
  
  while(reader.GetNextAlignment(al1)){
    delim_ind=whichDelim(al1.Name);
    
    if((i%1000000)==0) Rprintf("processing %d'th read\n",i);
    
    // Rprintf("processing %d'th read\n",i);

    /**
     * keepIt1 = 1 if the first read of the paired-end fragment is to be kept
     * keepIt2 = 1 if the second read of the paired-end fragment is to be kept 
     */
    keepIt1=1; 
    
    if(i++ < *skip){keepIt1=0;}
    
    if(al1.IsPaired() && al1.IsMateMapped()){
      // assume the reads are sorted by name, so the next read is the mate
      if(reader.GetNextAlignment(al2)){
        if((i%1000000)==0) Rprintf("processing %d million-th read\n",i/1000000);
        
        isSame=((split(al1.Name,DELIM_LOOKUP[delim_ind])[0]).compare(split(al2.Name,DELIM_LOOKUP[delim_ind])[0]) == 0);
        
        if (! isSame) {
          cout<<al1.Name<<endl;
          cout<<al2.Name<<endl;
          error("expect a paired-read, but see something else\n");
        }
        
        keepIt2=1;
        if(i++ < *skip){keepIt2=0;}
      }else{
        keepIt2=0;
      }
    }else{
      keepIt2=0;
    }
    
    if ((keepIt1==1) &&
        (((al1.AlignmentFlag & flag2kp) != flag2kp) || (al1.AlignmentFlag & flag2rm))){
      skip0 += 1;
      keepIt1 = 0;
    }
    
    if ((keepIt2==1) &&
        (((al2.AlignmentFlag & flag2kp) != flag2kp) || (al2.AlignmentFlag & flag2rm))){
      skip0 += 1;
      keepIt2 = 0;
    }
    
    if(keepIt1==1){
      avg=0;
      pBaseQualities = (char*)al1.Qualities.data();
      
      for(j=0;j<al1.Length;j++){
        avg+=pBaseQualities[j];
      }
      
      favg=(double)avg/(al1.Qualities).size()-phred;//replace with variable phred
      
      if(favg < avgQ){
        skip1 += 1;
        keepIt1=0;
      }
      
    }
    
    if(keepIt2==1){
      avg=0;
      pBaseQualities = (char*)al2.Qualities.data();
      
      for(j=0;j<al2.Length;j++){
        avg+=pBaseQualities[j];
      }
      
      favg=(double)avg/(al2.Qualities).size()-phred;//replace with variable phred
      
      if(favg<avgQ){
        skip1 += 1;
        keepIt2=0;
      }
    }
    
    if((keepIt1==0) && (keepIt2==0)){continue;}
    
    if((keepIt1==0) && (keepIt2==1)){
      keepIt1=1;
      keepIt2=0;
      al1=al2;
    }
    
    if((keepIt1!=0) && (keepIt2==0)){
      al1.MateRefID    = -1;
      al1.MatePosition = -1;
      al1.InsertSize   =  0;
      al1.SetIsMateMapped(false); //BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
      al1.SetIsPaired(false);     //BAM_ALIGNMENT_PAIRED=0x0001	  
    
      // we set offset = j + 1 because later we calculte the location by
      // al1.Position + offset1, and al1.Position with starting index 0
      for (j=0; j<al1.Length; j++) { offset1[j] = j+1; }
      
      /* seqPosj is the current position in the sequence */
      seqPosj=0;    
      cigarIter = al1.CigarData.begin();
      cigarEnd  = al1.CigarData.end();
      
      for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        const CigarOp& op = (*cigarIter);
        
        switch (op.Type) {
          case (Constants::BAM_CIGAR_MATCH_CHAR)    :
          case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
            seqPosj += op.Length;
            break;
          
          case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
            // because the postion of the sequence is the
            // 1-based leftmost mapping POSition of the first matching base.
            // need to hande the case where a mistmach happen at the beginning
            if (seqPosj==0) {
              for(j=(seqPosj+op.Length);j<al1.Length;j++){
                offset1[j]-= op.Length;
              }
            }
            seqPosj += op.Length;
            break;

          case (Constants::BAM_CIGAR_DEL_CHAR)      :
          case (Constants::BAM_CIGAR_REFSKIP_CHAR)  :
            // a sequence should not start with deletion or skipping
            // so we do not handle such cases here
            for (j=seqPosj; j<al1.Length; j++) {
              offset1[j] += op.Length; 
            }
            break;
            
          case (Constants::BAM_CIGAR_INS_CHAR)      :
          case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
            for (j=seqPosj; j<(seqPosj+op.Length); j++) {
              offset1[j] = -1;
            }
            for(j=(seqPosj+op.Length);j<al1.Length;j++){
              offset1[j]-= op.Length;
            }
            seqPosj  += op.Length;
            break;
        
          // hard clipping (clipped sequences NOT present in SEQ)
          // need to do nothing because the hard clipped sequence
          // is not in the reference, nor in the read sequence (SEQ)
          case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
            break;
        
          // shouldn't get here
          default:
            cout << al1.Name << endl;
            Rprintf("%c, %d\n", op.Type, op.Length);
            error("sorry, did not expect to see this character in cigar\n");
        }
      }
      
      k1 = k2 = k3 = 0;	
      sprintf(currentchr,"%s",chrs[al1.RefID].c_str());
      
      for(j=0; j<al1.Length; j++){
        /**
         * consider only the cases when we don't have an insertion
         * and quality of the read meats the criteria for a snp to be considered
         */
        if((offset1[j]!=-1) && ((int)al1.Qualities[j] >= (phred + snpQ))){
          currentpos=al1.Position+offset1[j];
          /* search the current read */
          sprintf(key, "%s.%ld", currentchr, currentpos);
          // Rprintf("%s %ld %ld\n", key, al1.Position, offset1[i]);
          
          if(hash_map->find(key, &alleles)) {
			       /**
             * make sure this snp is not methylated
             */
            nSNP += 1;
            if(toupper(al1.QueryBases[j]) == alleles[0]){
              pos1[k1++]=currentpos;
            }else if(toupper(al1.QueryBases[j]) == alleles[1]){
              pos2[k2++]=currentpos;
            }else {
              n_unrecognized += 1;
              pos3[k3++]=currentpos;
			        // Rprintf("observed allele %c, expect %c or %c\n", al1.QueryBases[j], alleles[0], alleles[1]);
            }
          }

        }else{
		      skip2 +=1;
		    }	
      }
      
      kall=k1+k2+k3;
      
      if(kall>0){
        prop1=((double)k1)/kall;
        prop2=((double)k2)/kall;
        
        if((prop1>propcut)||(prop2>propcut)){
          if(prop1>propcut){
            writer1.SaveAlignment(al1);
          }else{
            writer2.SaveAlignment(al1);
          }
        }else{
          // Rprintf("observed an inconsistent read: found snps of allele1 %d, allele2 %d & another allele %d, need to skip it.\n",k1,k2,k3);
          writer3.SaveAlignment(al1);
        }
      }
      
    }else{ // keepIt1!=0 and keepIt2!=0
      
      for (j=0;j<al1.Length;j++){offset1[j] = j+1;}
      for (j=0;j<al2.Length;j++){offset2[j] = j+1;}
      
      seqPosj=0;    
      cigarIter = al1.CigarData.begin();
      cigarEnd  = al1.CigarData.end();
      
      for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        const CigarOp& op = (*cigarIter);
        switch (op.Type) {
          case (Constants::BAM_CIGAR_MATCH_CHAR)    :
          case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
            seqPosj += op.Length;
            break;
          
          case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
            // because the postion of the sequence is the
            // 1-based leftmost mapping POSition of the first matching base.
            // need to hande the case where a mistmach happen at the beginning
            if (seqPosj==0) {
              for(j=(seqPosj+op.Length);j<al1.Length;j++){
                offset1[j]-= op.Length;
              }
            }
            seqPosj += op.Length;
            break;

          case (Constants::BAM_CIGAR_DEL_CHAR) :
          case (Constants::BAM_CIGAR_REFSKIP_CHAR)  :
            for (j=seqPosj; j<al1.Length; j++) { 
              offset1[j] += op.Length; 
            }
            break;
            
          case (Constants::BAM_CIGAR_INS_CHAR)      :
          case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
            for (j=seqPosj; j<(seqPosj+op.Length); j++) { 
              offset1[j] = -1;
            }
            for(j=(seqPosj+op.Length);j<al1.Length;j++){
              offset1[j]-= op.Length;
            }
            seqPosj  += op.Length;
            break;
        
          // hard clipping (clipped sequences NOT present in SEQ)
          // need to do nothing because the hard clipped sequence
          // is not in the reference, nor in the read sequence (SEQ)
          case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
            break;

          // shouldn't get here
          default:
            cout << al1.Name << endl;
            Rprintf("%c, %d", op.Type, op.Length);
            error("sorry, did not expect to see this character in cigar\n");
          }
	    }
      
      seqPosj=0;    
      cigarIter = al2.CigarData.begin();
      cigarEnd  = al2.CigarData.end();
      
      for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        const CigarOp& op = (*cigarIter);
        switch (op.Type) {
          case (Constants::BAM_CIGAR_MATCH_CHAR)    :
          case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
            seqPosj += op.Length;
            break;
        
          case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
            // because the postion of the sequence is the
            // 1-based leftmost mapping POSition of the first matching base.
            // need to hande the case where a mistmach happen at the beginning
            if (seqPosj==0) {
              for(j=(seqPosj+op.Length);j<al1.Length;j++){
                offset1[j]-= op.Length;
              }
            }
            seqPosj += op.Length;
            break;

          case (Constants::BAM_CIGAR_DEL_CHAR) :
          case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
            for (j=seqPosj; j<al2.Length; j++) { 
              offset2[j] += op.Length; 
            }
            break;
          case (Constants::BAM_CIGAR_INS_CHAR)      :
          case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
            for (j=seqPosj; j<(seqPosj+op.Length); j++) { 
              offset2[j] = -1;
            }
            for(j=(seqPosj+op.Length);j<al2.Length;j++){
              offset2[j]-= op.Length;
            }
            seqPosj  += op.Length;
            break;
            
          // hard clipping (clipped sequences NOT present in SEQ)
          // need to do nothing because the hard clipped sequence
          // is not in the reference, nor in the read sequence (SEQ)
          case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
            break;
            
          // shouldn't get here
          default:
            cout << al1.Name << endl;
            Rprintf("%c, %d", op.Type, op.Length);
            error("sorry, did not expect to see this character in cigar\n");
          }
      }
      
      k1 = k2 = k3 = 0;	
	    sprintf(currentchr,"%s",chrs[al1.RefID].c_str());
      
      for(j=0; j<al1.Length; j++){
        /**
         * consider only the cases when we don't have an insertion
		     * and quality of the read meats the criteria for a snp to be considered
         */
        
        if((offset1[j]!=-1)&&((int)al1.Qualities[j] >= (phred + snpQ))){
          currentpos=al1.Position+offset1[j];
          /* search the current read */
          sprintf(key, "%s.%ld", currentchr, currentpos);
          
          if(hash_map->find(key, &alleles)) {
            /**
             * make sure this snp is not methylated
             */
            nSNP += 1;
            if(toupper(al1.QueryBases[j]) == alleles[0]){
              pos1[k1++]=currentpos;
            }else if(toupper(al1.QueryBases[j]) == alleles[1]){
              pos2[k2++]=currentpos;
            }else {
              n_unrecognized += 1;
              pos3[k3++]=currentpos;

              // Rprintf("observed allele %c, expect %c or %c\n", al1.QueryBases[j], alleles[0], alleles[1]);
            }
          }
        }else{
          skip2 +=1;
        }					
      }
      
      for(j=0; j<al2.Length; j++){
        /**
         * consider only the cases when we don't have an insertion
		     * and quality of the read meats the criteria for a snp to be considered
         */
        
        if((offset2[j]!=-1)&&((int)al2.Qualities[j] >= (phred + snpQ))){
          currentpos=al2.Position+offset2[j];
          sprintf(key, "%s.%ld", currentchr, currentpos);
          
          if(hash_map->find(key, &alleles)){
			       /**
             * make sure this snp is not methylated
             */
            nSNP += 1;
            
            if(toupper(al2.QueryBases[j]) == alleles[0]){			  
              for(isAccounted=0;isAccounted<k1;isAccounted++){
                if(pos1[isAccounted]==currentpos)
                  offset2[j]=-1;
              }
              if(offset2[j]!=-1)pos1[k1++]=currentpos;
              
            }else if(toupper(al2.QueryBases[j]) == alleles[1]){
              for(isAccounted=0;isAccounted<k2;isAccounted++){
                if(pos2[isAccounted]==currentpos)
                  offset2[j]=-1;
              }
              if(offset2[j]!=-1)pos2[k2++]=currentpos;
              
            }else {
              n_unrecognized += 1;
              
              for(isAccounted=0;isAccounted<k3;isAccounted++){
                if(pos3[isAccounted]==currentpos)
                  offset2[j]=-1;
              }
              if(offset2[j]!=-1)pos3[k3++]=currentpos;
            
              // Rprintf("observed allele %c, expect %c or %c\n", al2.QueryBases[j], alleles[0], alleles[1]);
            }
          }
        }else{
          skip2 +=1;
        }	
      }
      
      kall=k1+k2+k3;
      
      if(kall>0){
        prop1=((double)k1)/kall;
        prop2=((double)k2)/kall;
        
        if((prop1>propcut)||(prop2>propcut)){
          if(prop1>propcut){
            writer1.SaveAlignment(al1);
            writer1.SaveAlignment(al2);
          }else{
            writer2.SaveAlignment(al1);
            writer2.SaveAlignment(al2);
          }
        }else{
          // Rprintf("observed an inconsistent read: found snps of allele1 %d, allele2 %d & another allele %d, need to skip it.\n",k1,k2,k3);
          writer3.SaveAlignment(al1);
          writer3.SaveAlignment(al2);
        }
      }
    }
    
  }
  
  reader.Close();
    
  writer1.Close();
  writer2.Close();
  writer3.Close();

  delete hash_map;
  free(key);
  free(alleles);
  
  Rprintf("total %d lines processed\n",i);
}
}// extern "C"
