my $finput  = "Homo_sapiens.GRCh37.66.gtf";
my $foutput = "Homo_sapiens.GRCh37.66.exon.gtf";


open(FIN, $finput) or die "Cannot find file $finput \n";
open(FOUT, ">".$foutput) or die "Cannot create file $foutput\n";

while(<FIN>) {
    if($_ =~ /^(\d|X|Y)/){
    	@values = split('\t', $_);
    	
    	if($values[2] eq "exon"){
    		print FOUT "chr".$_;
    	}
    }
}
close(FIN);
close(FOUT);
