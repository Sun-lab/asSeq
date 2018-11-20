/**********************************************************************
 *
 * getUnique_bam.cc
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2011, Vasyl Zhabotynsky, UNC-CH
 *
 * first written Jul 27, 2011
 * 
 * modified by Chong Jin 2016-06-28
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 **********************************************************************/

//#include <stdlib.h>
//#include <math.h>
//#include <time.h>
#include <R.h>
#include "BamMultiReader.h"
#include "BamWriter.h"
#include "bamtools_sort.h"
using namespace BamTools;
#include <vector>
#include <string>
using namespace std;
#include "alName.h"

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

void getUnique_bam(char **Rinput, int* Rnum_inp, char **Routput, int *fixmate)
{
  int delim_ind=0;
  int num_inp  = *Rnum_inp;
  int i, flag, chri=0;

  string output;
  output=string(Routput[0]);
 
  BamMultiReader reader;
  BamWriter writer;

  const char* strand="+-";
  //strand[0]='+';
  //strand[1]='-';
	
  std::vector<std::string> input,chrs;
  std::string temp;
  
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
  
  if(!writer.Open(output,header,references)){
	  Rprintf("Could not open %s output file.\n",output.c_str());
	  return;
  }

  SamSequenceConstIterator begin = header.Sequences.ConstBegin();
  SamSequenceConstIterator iter  = begin;
  SamSequenceConstIterator end   = header.Sequences.ConstEnd();
  
  for ( ; iter != end; ++iter ) {
    const SamSequence& currentSeq = (*iter);
    chrs.push_back(currentSeq.Name);
    chri++;
  }

  int fixedPair   = 0;
  int skipPair  = 0;
  int skipOne   = 0;
  int skipMulti = 0;
  int writPair  = 0;
  int writOne   = 0;
  
  BamAlignment al1,al2,tmp;
  
  //al1.Name="";
  //cout<<al1.Name<<endl;
  int notNullAlignment=0,isSame=0;
  int numReads=1;
  
  notNullAlignment=reader.GetNextAlignment(al1);
  
  i=1;
  delim_ind=whichDelim(al1.Name);
  cout<<"delim: "<<DELIM_LOOKUP[delim_ind]<<endl;

  //cout<<al1.Name<<" "<<numReads<<" "<<delim_ind<<" "<<DELIM_LOOKUP[delim_ind]<<endl;
  while(notNullAlignment){        
    //if(i>100)break;
    // Find all adjacent alignments with same qname.
    do{
      notNullAlignment=reader.GetNextAlignment(tmp);

      //cout<<tmp.Name<<" "<<numReads<<endl;
      //read in alignments,
      //if 1 or 2 process
      //otherwise make i=1, and reread alignments
      if(notNullAlignment){
        i++;
        if((i%1000000)==0)Rprintf(".");
        if((i%10000000)==0)Rprintf("\nprocessing %d million-th line\n", i/1000000);

        //cout<<split(al1.Name,DELIM_LOOKUP[delim_ind])[0]<<" "<<split(tmp.Name,DELIM_LOOKUP[delim_ind])[0]<<endl;
        isSame=((split(al1.Name,DELIM_LOOKUP[delim_ind])[0]).compare(split(tmp.Name,DELIM_LOOKUP[delim_ind])[0]) == 0);
        if(isSame){
          al2=tmp;
          numReads++;		  
        }
        //cout<<al1.Name<<" "<<tmp.Name<<" "<<isSame<<" "<<numReads<<endl;
      }else{
          break;
      }      
      //cout << numReads<<endl;
    }while(isSame);
    
    if(numReads==1){
      al1.MateRefID    = -1;
      al1.MatePosition = -1;
      al1.InsertSize   = 0;
      al1.SetIsMateMapped(true);//~BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
      al1.SetIsPaired(false);//BAM_ALIGNMENT_PAIRED=0x0001	 
      writer.SaveAlignment(al1);
      writOne++;
    }else if(numReads==2){
      // samtools fixmate for qname-sorted reads. -Chong
      if (*fixmate && (!al1.IsPaired() || !al2.IsPaired())) {
        al1.SetIsPaired(true);
        al2.SetIsPaired(true);
        al1.MateRefID = al2.RefID; // fix RNEXT
        al2.MateRefID = al1.RefID;
        al1.MatePosition = al2.Position; // fix PNEXT
        al2.MatePosition = al1.Position;
        fixedPair+=numReads;
      }
      if((al1.IsPaired())&&(al2.IsPaired())){
        if((al1.IsMapped())&&(al2.IsMapped())){
          al1.SetIsMateMapped(true);
          al2.SetIsMateMapped(true);
          writer.SaveAlignment(al1);
          writer.SaveAlignment(al2);
          writPair++;
        }else if((al1.IsMapped())&&(!al2.IsMapped())){
          al1.MateRefID    = -1;
          al1.MatePosition = -1;
          al1.InsertSize   = 0;
          // "If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80." -- Sam format specification
          // Unset 0x8 to avoid "Mate unmapped flag should not be set for unpaired reads" error in picard.
          al1.SetIsMateMapped(true);//~BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
          al1.SetIsPaired(false);//BAM_ALIGNMENT_PAIRED=0x0001	 
          writer.SaveAlignment(al1);		  
          skipOne++;
          writOne++;
        }else if((!al1.IsMapped())&&(al2.IsMapped())){
          al2.MateRefID    = -1;
          al2.MatePosition = -1;
          al2.InsertSize   = 0;
          // Unset 0x8 to avoid "Mate unmapped flag should not be set for unpaired reads" error in picard.
          al2.SetIsMateMapped(true);//~BAM_ALIGNMENT_MATE_UNMAPPED=0x0008
          al2.SetIsPaired(false);//BAM_ALIGNMENT_PAIRED=0x0001	 
          writer.SaveAlignment(al2);		  
          skipOne++;
          writOne++;
        }else{
          skipPair+=numReads;
        }
      }else{
        skipPair+=numReads;
      }
    }else{
      skipMulti+=numReads;
    }
    
    al1=tmp;
    numReads=1;
    //cout << numReads<<endl;
  }
  
  Rprintf("%d paired reads with same qnames have isPaired flag fixed\n", fixedPair);
  Rprintf("keep %d reads that are not paired\n", writOne);
  Rprintf("keep %d reads that are paired\n", writPair*2);
  Rprintf("%d single reads are skipped because they are not mapped\n", skipOne);
  Rprintf("%d paired reads are skipped because they are not mapped\n", skipPair);
  Rprintf("%d reads that appear more than twice are skipped\n", skipMulti);
  Rprintf("total %d lines processed\n", i);


  reader.Close();
  writer.Close();
}
}// extern "C"
