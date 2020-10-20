#!/usr/bin/python
# -*- coding: utf-8 -*-
# 11 fevrier 2019
# Compare les reads 1 pour identifier d'eventuels 'twins' (sites apparemment differents mais qui correspondraient a un seul et meme site ; multiplication liee a des erreurs de mapping sur des sequences similaires)
# Par rapport a twin.py : compare les sisters (evalues par soniclength) et non les shear sites.

import sys #twin.py [fichier d'entrée : count_SFV11c_removed.txt] [fichier de sortie : twins_SFV11c.txt] [genome]
import os
import os.path
from os import path

file = open(sys.argv[1],'r')
sites = file.readlines()
file.close()

CLUS = []
inclus = []

## Recherche de twins eventuels ##
for i in range(1,len(sites)-1) :
	li = sites[i].strip().split()
	
	file = open('seqi.fa','w')
	file.write('>'+li[1]+'_'+li[2]+'\n'+li[8]+'\n')
	#file.write('>'+li[1]+'_'+li[2]+'\n'+li[7]+'\n')
	file.close()
	
	#print 'i :',li[1],' ',li[2]
	
	group = [li]
	for j in range(i+1,len(sites)) :
		lj = sites[j].strip().split()
		
		file = open('seqj.fa','w')
		file.write('>'+lj[1]+'_'+lj[2]+'\n'+lj[8]+'\n')
		#file.write('>'+lj[1]+'_'+lj[2]+'\n'+lj[7]+'\n')
		file.close()
		
		#print 'j :',lj[1],' ',lj[2]
		
		os.popen("blastn -query seqi.fa -subject seqj.fa -dust no -out bl2seq.out -outfmt 6")
		
		if path.exists("bl2seq.out") :
			file = open('bl2seq.out','r')
			out = file.readlines()
			file.close()
			
			if len(out) != 0 :
				lo = out[0].split()
				if float(lo[2]) >= 95 and int(lo[3]) >= 25 :
					#print out[0]
					group.append(lj)
	
	if len(group) > 1 :
		test = 0
		for g in group :
			if g[0]+'_'+g[1] in inclus :
				test += 1
			else :
				inclus.append(g[0]+'_'+g[1])
		
		if test != len(group) :
			#print group
			CLUS.append(group)

print CLUS

## Pour chaque cluster, blast des sequences contre le genome et recup du meilleur hit ##

file = open('seq.fa','w')

ecrit = []

for c in CLUS :
	for i in c :
		if i[1]+'_'+i[2] not in ecrit :
			ecrit.append(i[1]+'_'+i[2])
			file.write('>'+i[1]+'_'+i[2]+'\n'+i[8]+'\n')
			#file.write('>'+i[1]+'_'+i[2]+'\n'+i[7]+'\n')

file.close()

os.popen('blastn -query seq.fa -subject '+sys.argv[3]+' -dust no -out blastn.out -outfmt 6')

file = open('blastn.out','r')
res = file.readlines()
file.close()

RES = {}
for r in res :
	lr = r.strip().split()
	site = lr[0]
	if site not in RES :
		RES[site] = r

file = open(sys.argv[2],'w')
n = 1
for c in CLUS :
	for i in c :
		print i
		if i[1]+'_'+i[2] in RES :
			file.write('clus'+str(n)+'\t'+RES[i[1]+'_'+i[2]])
	n += 1
file.close()

## Modif du fichier de count ##
# Inversion du dico des clusters
SIC = {}
n = 1
CLUSOUT = []

for c in CLUS :
	max_shear = 0
	max_sister = 0.0
	max_perc = 0
	max_read = 0
	seqmax = ''
	totread = 0
	all_sites = ''
	
	for i in c :
		all_sites = all_sites+i[1]+'_'+i[2]+'_'+i[3]+'_'+i[4]+'_'+i[5]+'_'+i[6]+' '+i[7]+' '
		
		if i[1]+'_'+i[2] not in SIC :
			SIC[i[1]+'_'+i[2]] = 'clus'+str(n)
		else :
			SIC[i[1]+'_'+i[2]] = SIC[i[1]+'_'+i[2]] + '_clus'+str(n)
		
		totread += int(i[7])
		
		#if int(i[4]) > max_shear :
		if float(i[5]) > max_sister :
			max_shear = int(i[4])
			max_sister = float(i[5])
			max_perc = i[6]
			max_read = int(i[7])
			seqmax = i[8]
		#elif int(i[4]) == max_shear :
		elif float(i[5]) == max_sister :
			if int(i[7]) > max_read :
				max_shear = int(i[4])
				max_read = int(i[7])
				max_perc = i[6]
				seqmax = i[8]
		
	
	CLUSOUT.append(i[0]+'\tcluster'+str(n)+'\tNA\tNA\t'+str(max_shear)+'\t'+str(max_sister)+'\t'+max_perc+'\t'+str(totread)+'\t'+seqmax+'\t'+all_sites)
		
	n += 1
	
print sites

file = open(sys.argv[1].replace('.txt','_inclust.txt'),'w')
file.write(sites[0].strip()+'\tClusters\n')

for i in range(1,len(sites)) :
	li = sites[i].strip().split('\t')
	file.write(sites[i].strip()+'\t')
	if li[1]+'_'+li[2] in SIC :
		file.write(SIC[li[1]+'_'+li[2]]+'\n')
	else :
		file.write('no_twin\n')

file.close()

# Recalcul du nombre total de shear sites et  sisters #
#totshear = 0 
totsister = 0.0

for i in range(1,len(sites)) :
	li = sites[i].strip().split('\t')
	if li[1]+'_'+li[2] not in SIC :
		#totshear += int(li[4])
		totsister += float(li[5])

for c in CLUSOUT :
	#totshear += int(c.split('\t')[4])
	totsister += float(c.split('\t')[5])
		

file = open(sys.argv[1].replace('.txt','_clustered.txt'),'w')
file.write(sites[0].strip()+'\tSites_in_clust\n')

for i in range(1,len(sites)) :
	li = sites[i].strip().split('\t')
	if li[1]+'_'+li[2] not in SIC :
		file.write(li[0]+'\t'+li[1]+'\t'+li[2]+'\t'+li[3]+'\t'+li[4]+'\t'+li[5]+'\t'+str(float(li[5])*100/float(totsister))+'\t'+li[7]+'\t'+li[8]+'\tno_twin\n')
		
for c in CLUSOUT :
	lc = c.split('\t')
	file.write(lc[0]+'\t'+lc[1]+'\t'+lc[2]+'\t'+lc[3]+'\t'+lc[4]+'\t'+li[5]+'\t'+str(float(lc[5])*100/float(totsister))+'\t'+lc[7]+'\t'+lc[8]+'\t'+lc[9]+'\n')

file.close()



