#!/usr/bin/python
# -*- coding: utf-8 -*-


import os

results = "results.blast"
ref = "meiosis-related-genes.fa"


ld = {}
fasta = "/media/antoine/DATA3/vagawild/final_01/meiosis/meiosis.fa"
f = open(fasta, 'r')
for line in f :
	if line == "" :
		continue
	if line[0] == ">" :
		ld[line.strip()[1:]] = 0
	else :
		ld[ list(ld.keys())[-1] ] += len(line.strip())
f.close()
print(ld)
print('')


class CTG :
	def __init__(self) :
		self.genes = {}

	def addGene(self, gene, new_region) :
		if gene not in self.genes.keys() :
				self.genes[gene] = [new_region]
		else :
			overlap = False
			for stored_region in self.genes[gene] :
				if stored_region.isOverlapping(new_region) :
					overlap = True					
					break
			if not overlap :
				self.genes[gene].append(new_region)

	def printhits(self) :
		for gene in self.genes :
			toprint = gene + ": "
			for region in self.genes[gene] :
				toprint += str(region.start) + "-" + str(region.end) + "  "
			print(toprint)

class REG :
	def __init__(self, s, e) :
		self.start = s
		self.end = e
		if self.end < self.start :
			self.start, self.end = self.end, self.start

	def isOverlapping(self, region) :
		if region.start <= self.start and region.end >= self.end :
			return True
		elif region.end > self.start and region.end < self.end :
			return True
		elif region.start > self.start and region.start < self.end :
			return True
		else :
			return False
			


contigs = {}

f = open(results, 'r')
for line in f :
	line = line.strip()
	s = line.split('\t')
	hit = s[1]
	length = int(s[3])
	# CHANGE HIT IN DICT IF BETTER THAN PREVIOUS
	pident = float(s[2])
	bitscore = float(s[-1])
	evalue = float(s[7])
	completeness = float(length/ld[hit])
	start = int(s[5])
	end = int(s[6])
	chromosome = s[0]

	if end < start :
		start, end = end, start

	if chromosome not in contigs.keys() :
		contigs[chromosome] = CTG()
		if evalue < 0.000001 and bitscore > 50.0 :
			region = REG(start, end)
			gene = hit
			contigs[chromosome].addGene(gene, region)
	else :
		if evalue < 0.000001 and bitscore > 50.0 :
			region = REG(start, end)
			gene = hit
			contigs[chromosome].addGene(gene, region)


hits_count = {}
for ctg in contigs :
	print(ctg)
	contigs[ctg].printhits()
	for gene in contigs[ctg].genes :
		try :
			hits_count[gene] += 1
		except :
			hits_count[gene] = 1
	print("")

for gene in hits_count :
	print(gene, hits_count[gene])

f.close()
