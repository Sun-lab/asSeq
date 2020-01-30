/*
 * Parse a binary 1000 genomes vcf file into two-column text output like:
 * variant_number\tsample_number\n
 * variant_number\tsample_number\n
 *  ...          \t ...
 * that is, it only shows whether or not a variant occured not what it is.
 * NOTE!! This program ignores phasing and assumes specific 1000 genomes vcf format.
 * Example use:
 * cc -O2 parse.c
 * zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz  | cut  -f "10-" | ./a.out 
 * updated to add phasing information
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv)
{
  int k = 0, t, v;
  char *entry;
  char *line = NULL;
  size_t size, n;
  while((n = getline(&line, &size, stdin)) != -1)
  {
    if(line[0] == '#') continue;
    k++;
    t = 1;
    entry = strtok(line, "\t");
    while(entry != NULL)
    {
      if(strncmp("0|0", entry, 3) != 0){
        v=1;
         /*AA - 0, AB - 1, BA - 2, BB - 3*/
        if(strncmp(".|.", entry, 3) == 0)v=-1;
        if(strncmp("0|1", entry, 3) == 0)v=1;
        if(strncmp("1|0", entry, 3) == 0)v=2;
        if(strncmp("1|1", entry, 3) == 0)v=3;
        printf("%d\t%d\t%d\n", k, t, v);      
      }
      entry = strtok(NULL, "\t");
      t++;
    }
  }
  if(line) free(line);
  return 0;
}