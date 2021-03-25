/**********************************************************************
 *
 * asCounts_sam.c
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2011, Vasyl Zhabotynsky, UNC-CH
 *
 * first written Jul 27, 2011
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 **********************************************************************/

#include "BamMultiReader.h"
#include "BamWriter.h"
#include "bamtools_sort.h"
using namespace BamTools;

extern "C" {
void sort_bam(char **Rparams, int* Rnum_inp)
{  
  int num_inp;
  num_inp=*Rnum_inp;

  SortTool tool;
  tool.Run(num_inp, Rparams);
}
}// extern "C"
