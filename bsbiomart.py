import bioservices as bs
import sys

def gethumanshortname(outfile):
	#Get all human ENSEMBL / short name relationships
	s = bs.BioMart(host = 'www.ensembl.org')
	s.add_dataset_to_xml('hsapiens_gene_ensembl')
	s.add_attribute_to_xml('ensembl_gene_id')
	s.add_attribute_to_xml('wikigene_name')
	xml = s.get_xml()
	res = s.query(xml)

	res = res.split('\n')
	d = {} #{ENSG : shortname}
	for pair in res:
		pair = pair.split('\t')
		if len(pair) == 2:
			#If there is a shortname
			ensembl = str(pair[0])
			shortname = str(pair[1])
			d[ensembl] = shortname
		elif len(pair) == 1:
			#If there is not a shortname
			ensembl = str(pair[0])
			d[ensembl] = None

	with open(outfile, 'w') as f:
		for gene in d:
			if d[gene]:
				f.write(gene + '\t' + d[gene] + '\n')


def mousetohuman(genelist, outfile):
	#Given a list of ENSEMBL mouse gene IDs, output file of human orthologs

	#First get all mouse/human ortholog relationships
	s = bs.BioMart(host = 'www.ensembl.org')
	s.add_dataset_to_xml('mmusculus_gene_ensembl')
	s.add_attribute_to_xml('ensembl_gene_id')
	s.add_attribute_to_xml('hsapiens_homolog_ensembl_gene')
	xml = s.get_xml()
	res = s.query(xml)

	res = res.split('\n')
	d = {} #{ENSMUSG : ENSG}
	for pair in res:
		pair = pair.split('\t')
		if len(pair) == 2:
			#If there is a human ortholog
			mousegene = str(pair[0])
			humangene = str(pair[1])
			d[mousegene] = humangene
		elif len(pair) == 1:
			#If there is not a human homolog
			mousegene = str(pair[0])
			d[mousegene] = None

	with open(outfile, 'w') as f, open(genelist, 'r') as i:
		genes = []
		for line in i:
			line = line.strip()
			genes.append(line)
		for gene in d:
			if d[gene] and gene in genes:
				f.write(d[gene] + '\n')



if __name__ == '__main__':
	mousetohuman(sys.argv[1], sys.argv[2])
