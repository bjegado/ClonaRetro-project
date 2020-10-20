#!/usr/bin/python
# -*- coding: utf-8 -*-
# Par Magali Naville, IGFL Lyon
# 4 mai 2018
# Elimine les sites dupliqués entre individus. Les sites dupliqués sont attribues à l'espèce présentant le plus grand nombre de sites de coupure ou de duplicats PCR associes, pour un nombre de sites de coupure > 6.

## Doublons ##
#file = open('/chemin/doublons_rmart.list','r')
file = open('doublons_rmart.list','r')
doub = file.readlines()
file.close()

DUP = {}
IND = {}
for d in doub :
	ld = d.split('\t')
	DUP[ld[0]] = {} # dico, pour chaque site, des indivdus avec le nombre de sites de coupures et le nombre de duplicats PCR
	IND[ld[0]] = '' # liste des individus presentant le site
	for i in range(1,len(ld)-1) :
		IND[ld[0]] = IND[ld[0]] + ld[i].split()[0] + ' '
		DUP[ld[0]][ld[i].split()[0]] = [ld[i].split()[1], ld[i].split()[2]]
#print DUP

RETAINED = {} # Dico des individus retenus pour chaque site retenu
for i in IND :
	if IND[i] == 'SFV13a SFV13b ' or IND[i] == 'SFV11c SFV11t ' :
		RETAINED[i] = IND[i]
	else :
		best = ''
		SS = 6 # nombre minimal de sites de coupure requis pour garder le site
		reads = 0
		listind = IND[i].split()
		#print i
		#print listind
		#print DUP[i]
		for l in listind :
			#print l
			#print DUP[i][l]
			if int(DUP[i][l][0]) > SS :
				best = l
				SS = int(DUP[i][l][0])
			elif int(DUP[i][l][0]) == SS :
				if int(DUP[i][l][1]) > reads :
					best = l
					reads = int(DUP[i][l][1])
		#print best
		if best != '' :
			RETAINED[i] = best

print RETAINED


# Indiquer les fichiers pour chaque échantillon

barcode = ['3_SFV_co','4_SFV_co','10_SFV_co','12_SFV_co','14_SFV_co','20_SFV_co','21_SFV_co','24_SFV_co','25_SFV_co','26_SFV_co','27_SFV_co','29_SFV_co','32_SFV_co','33_SFV_co','36_SFV_co','38_SFV_co','41_SFV_co','42_SFV_co','44_SFV_co','JKT_SFV'] 

#['1_SFV_mono','3_SFV_mono','4_SFV_mono','5_SFV_mono','6_SFV_mono','7_SFV_mono','8_SFV_mono','10_SFV_mono','11_SFV_mono','12_SFV_mono','13_SFV_mono','14_SFV_mono','15_SFV_mono','16_SFV_mono','17_SFV_mono','18_SFV_mono','19_SFV_mono','3_SFV_co','4_SFV_co','10_SFV_co','12_SFV_co','14_SFV_co','20_SFV_co','21_SFV_co','24_SFV_co','25_SFV_co','26_SFV_co','27_SFV_co','29_SFV_co','32_SFV_co','33_SFV_co','36_SFV_co','38_SFV_co','41_SFV_co','42_SFV_co','44_SFV_co','JKT_SFV'] 
#['3_STLV1_co','4_STLV1_co','10_STLV1_co','12_STLV1_co','14_STLV1_co','20_STLV1_co','21_STLV1_co','24_STLV1_co','25_STLV1_co','26_STLV1_co','27_STLV1_co','29_STLV1_co','32_STLV1_co','33_STLV1_co','36_STLV1_co','38_STLV1_co','41_STLV1_co','42_STLV1_co','44_STLV1_co','STLV1_JKT']


for b in barcode :
	tot = 0 # nombre de sites éliminés
	TOT = 0 # nombre de sites restant
	RET = [] # Sites retenus
	shs = 0 # nombre total de sites de coupure
	
	#file = open('count_'+b+'_UMI_ShearSite.txt','r')
	file = open('R12_'+b+'_rmart.txt','r')
	site = file.readlines()
	file.close()
		
	for i in range(1,len(site)) :
		li = site[i].strip().split()
		ins = li[1]+' '+li[2]
		
		if ins not in DUP :
			RET.append(site[i])
			TOT += 1
			shs += int(li[6])
		elif ins in RETAINED :
			if b in RETAINED[ins] :
				RET.append(site[i])
				TOT += 1
				shs += int(li[6])
			else : 
				tot += 1
		else :
			tot += 1
	
	file = open('count_'+b+'_rmart_no-double.txt','w')
	file.write('Individu\tChrom\tInsertion_Site\tStrand\tShear_sites\tPerc_of_total_ShearSites\tTotal_duplicates(reads)\tUMIs\tSeq_Read1\n')
	for r in RET :
		lr = r.split('\t')
		file.write(lr[0]+'\t'+lr[1]+'\t'+lr[2]+'\t'+lr[3]+'\t'+lr[6]+'\t'+str(float(lr[6])*100/float(shs))+'\t'+lr[7]+'\t'+lr[4]+'\t'+lr[8])
	file.close()
	
	print b, ' gardés: ', TOT, ' éliminés: ', tot
	

