my $finput  = "Mus_musculus.NCBIM37.67.gtf";
my $foutput = "Mus_musculus.NCBIM37.67.updated.gtf";


open(FIN, $finput) or die "Cannot find file $finput \n";
open(FOUT, ">".$foutput) or die "Cannot create file $foutput\n";

while(<FIN>) {
  s/^MT/M/g;
  
	if($_ =~ /^(\d|X|Y|M)/){
		print FOUT "chr".$_;
	}
}
close(FIN);
close(FOUT);
