#!/usr/bin/python
# -*- coding: utf-8 -*-
# Par Magali Naville, IGFL, Lyon
# 16 janvier 2019 
# Vérifie les sites identifiés en reblastant les sequences adjacentes dans le génome pour vérifier qu'elles ne ressemblent pas a la fin du LTR.

import sys
# delete_artefact_sites.py count_SFV163.txt SFV163 /home/mnaville/Documents/Clonaretro/Genome_babouin/Papio_anubis.Panu_3.0.dna.toplevel.fa /home/mnaville/Documents/Clonaretro/Genome_babouin/SFV_fin_LTR.fa
import os
import subprocess

# Conversion des sites en fichier bed += 100 nts
print('Conversion des sites en fichier bed += 100 nts')
os.popen("nawk 'NR>1 {print $2\"\t\"$3-100\"\t\"$3+100}' "+sys.argv[1]+" > sites_"+sys.argv[2]+".bed")
#subprocess.Popen("nawk 'NR>1 {print $2\"\t\"$3-100\"\t\"$3+100}' "+sys.argv[1]+" > sites_"+sys.argv[2]+".bed").wait()

# Modif du Bed pour éviter les coordonnées négatives
file = open("sites_"+sys.argv[2]+".bed",'r')
pos = file.readlines()
file.close()

file = open("sites_"+sys.argv[2]+"_cor.bed",'w')
for p in pos :
	print(p)
	lp = p.strip().split()
	if lp[1][0] == '-' :
		file.write(lp[0]+'\t0\t'+lp[2]+'\n')
	else :
		file.write(p)

# Récupération des séquences dans le génome
print('Recup des sequences dans le genome')
subprocess.call(["python","/home/brice/Documents/Script_files_LMPCR/get_seq_bed.py",sys.argv[3],"sites_"+sys.argv[2]+"_cor.bed","sites_"+sys.argv[2]+".fa"])
#os.popen("python /home/brice/Documents/Script_files_LMPCR/get_seq_bed.py "+sys.argv[3]+" sites_"+sys.argv[2]+"_cor.bed sites_"+sys.argv[2]+".fa")

# Remplacement des espaces par des "_" (pour garder la position génomique entière lors du Blast)
print('Modif des noms de sequences')
os.popen("sed -i 's/ /_/g' sites_"+sys.argv[2]+".fa")

# Blast du LTR contre les séquences génomiques
print('Blast du LTR contre les sequences genomiques')
os.popen("blastn -query "+sys.argv[4]+" -subject sites_"+sys.argv[2]+".fa -out LTR_tab.out -evalue 3 -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -outfmt 6")
os.popen("blastn -query "+sys.argv[4]+" -subject sites_"+sys.argv[2]+".fa -out LTR_aln.out -evalue 3 -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2")

# Parsing de la sortie Blast : sélection des hits se terminant à la position 99 (al +/+) ou 100 (al +/-) du fragment génomique
print('Parsing de la sortie Blast pour lister les sites artefactuels')

file = open('LTR_tab.out','r')
hits = file.readlines()
file.close()

DEL = []

for h in hits :
	if h[0] != '#' :
		lh = h.strip().split()
		if lh[9] == '99' or lh[9] == '100' :
			DEL.append(lh[1].strip('_'))

#print DEL

# Suppression des sites dans le fichier count
print('Suppression des sites dans le fichier count')
file = open(sys.argv[1],'r')
sites = file.readlines()
file.close()

file = open(sys.argv[1].replace('.sites.txt','_rmart.txt'),'w')
file.write(sites[0])
for s in range(1,len(sites)) :
	ls = sites[s].strip().split()
	if ls[1]+'_'+str(int(ls[2])-100)+'_'+str(int(ls[2])+100) not in DEL :
		file.write(sites[s])

file.close()

# Listing des sites supprimes dans un fichier
file = open(sys.argv[1].replace('.sites.txt','_removed.txt'),'w')
for d in DEL :
	pos = d.split('_')
	file.write(pos[0]+'\t'+str(int(pos[1])+100)+'\n')
file.close()
	

print('Nombre de sites identifiés : ')
print(len(sites)-1)
print('Nombre de sites supprimés : ')
print(len(DEL))
print(str(len(DEL)*100/(len(sites)-1))+' %')
