echo "cut -f 3 Homo_sapiens.GRCh37.66.gtf | sort | uniq -c"
echo ""
cut -f 3 Homo_sapiens.GRCh37.66.gtf | sort | uniq -c 
echo ""

echo "cut -f 1 Homo_sapiens.GRCh37.66.gtf | sort | uniq -c"
echo ""
cut -f 1 Homo_sapiens.GRCh37.66.gtf | sort | uniq -c
echo ""

echo "cut -f 1 Homo_sapiens.GRCh37.66.exon.gtf | sort | uniq -c"
echo ""
cut -f 1 Homo_sapiens.GRCh37.66.exon.gtf | sort | uniq -c 
echo ""
