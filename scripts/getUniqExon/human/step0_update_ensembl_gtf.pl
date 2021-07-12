my $finput  = "Homo_sapiens.GRCh37.66.gtf";
my $foutput = "Homo_sapiens.GRCh37.66.updated.gtf";


open(FIN, $finput) or die "Cannot find file $finput \n";
open(FOUT, ">".$foutput) or die "Cannot create file $foutput\n";

while(<FIN>) {	
    if($_ =~ /^(\d|X|Y)/){
    	print FOUT "chr".$_;
    }
}
close(FIN);
close(FOUT);

