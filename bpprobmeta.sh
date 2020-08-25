
echo ROupUTR3
python BpProbMeta.py --fasta ../rouptranscripts.UTR3.fa --windowsize 100 --slidesize 10 --numberofbins 50 --region UTR3 --seqclass roup --outputfile BpprobMeta.txt

echo ROdownUTR3
python BpProbMeta.py --fasta ../rodowntranscripts.UTR3.fa --windowsize 100 --slidesize 10 --numberofbins 50 --region UTR3 --seqclass rodown --outputfile BpprobMeta.txt

echo ROunchangedUTR3
python BpProbMeta.py --fasta ../rounchangedtranscripts.UTR3.fa --windowsize 100 --slidesize 10 --numberofbins 50 --region UTR3 --seqclass rounchanged --outputfile BpprobMeta.txt

echo ROupUTR5
python BpProbMeta.py --fasta ../rouptranscripts.UTR5.fa --windowsize 100 --slidesize 10 --numberofbins 20 --region UTR5 --seqclass roup --outputfile BpprobMeta.txt

echo ROdownUTR5
python BpProbMeta.py --fasta ../rodowntranscripts.UTR5.fa --windowsize 100 --slidesize 10 --numberofbins 20 --region UTR5 --seqclass rodown --outputfile BpprobMeta.txt

echo ROunchangedUTR5
python BpProbMeta.py --fasta ../rounchangedtranscripts.UTR5.fa --windowsize 100 --slidesize 10 --numberofbins 20 --region UTR5 --seqclass rounchanged --outputfile BpprobMeta.txt

echo ROupCDS
python BpProbMeta.py --fasta ../rouptranscripts.CDS.fa --windowsize 100 --slidesize 10 --numberofbins 100 --region CDS --seqclass roup --outputfile BpprobMeta.txt

echo ROdownCDS
python BpProbMeta.py --fasta ../rodowntranscripts.CDS.fa --windowsize 100 --slidesize 10 --numberofbins 100 --region CDS --seqclass rodown --outputfile BpprobMeta.txt

echo ROunchangedCDS
python BpProbMeta.py --fasta ../rounchangedtranscripts.CDS.fa --windowsize 100 --slidesize 10 --numberofbins 100 --region CDS --seqclass rounchanged --outputfile BpprobMeta.txt