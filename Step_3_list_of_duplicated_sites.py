#!/usr/bin/python
# -*- coding: utf-8 -*-
# Par Magali Naville, IGFL Lyon
# 10 octobre 2017 
# Recence les eventuels sites 'doublons' entre individus (= erreur probable dans le code-barre de l'adaptateur).

# Indiquer les fichiers pour chaque échantillon
barcode = ['3_STLV1_co','4_STLV1_co','10_STLV1_co','12_STLV1_co','14_STLV1_co','20_STLV1_co','21_STLV1_co','24_STLV1_co','25_STLV1_co','26_STLV1_co','27_STLV1_co','29_STLV1_co','32_STLV1_co','33_STLV1_co','36_STLV1_co','38_STLV1_co','41_STLV1_co','42_STLV1_co','44_STLV1_co','STLV1_JKT']

#['1_SFV_mono','3_SFV_mono','4_SFV_mono','5_SFV_mono','6_SFV_mono','7_SFV_mono','8_SFV_mono','10_SFV_mono','11_SFV_mono','12_SFV_mono','13_SFV_mono','14_SFV_mono','15_SFV_mono','16_SFV_mono','17_SFV_mono','18_SFV_mono','19_SFV_mono','3_SFV_co','4_SFV_co','10_SFV_co','12_SFV_co','14_SFV_co','20_SFV_co','21_SFV_co','24_SFV_co','25_SFV_co','26_SFV_co','27_SFV_co','29_SFV_co','32_SFV_co','33_SFV_co','36_SFV_co','38_SFV_co','41_SFV_co','42_SFV_co','44_SFV_co','JKT_SFV'] 
#['3_STLV1_co','4_STLV1_co','10_STLV1_co','12_STLV1_co','14_STLV1_co','20_STLV1_co','21_STLV1_co','24_STLV1_co','25_STLV1_co','26_STLV1_co','27_STLV1_co','29_STLV1_co','32_STLV1_co','33_STLV1_co','36_STLV1_co','38_STLV1_co','41_STLV1_co','42_STLV1_co','44_STLV1_co','STLV1_JKT']


SITE = {}

for b in barcode :
	file = open('R12_'+b+'_rmart.txt','r')
	eff = file.readlines()
	file.close()
	
	for e in range(1,len(eff)) :
		le = eff[e].strip().split()
		if len(le) > 4 :
			if le[1]+' '+le[2] not in SITE :
				#SITE[le[1]+' '+le[2]] = [b+' '+le[3]+' '+le[5]]
				SITE[le[1]+' '+le[2]] = [b+' '+le[6]+' '+le[4]]
			else :
				#SITE[le[1]+' '+le[2]].append(b+' '+le[3]+' '+le[5])
				SITE[le[1]+' '+le[2]].append(b+' '+le[6]+' '+le[4])

#file = open('doublons.list','w')
file = open('/home/brice/Bureau/doublons_rmart.list','w')
for s in SITE :
	if len(SITE[s]) > 1 :
		file.write(s + '\t')
		for i in SITE[s] :
			file.write(i + '\t')
		file.write('\n')
file.close()

print (len(SITE))

