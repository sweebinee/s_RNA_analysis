ls -d ERR* > N_171_sample.txt

import os, re

normal='N_171'
sample_pattern='ERR315'

maindir='/storage2/Project/TCGA_fusion/STARFusion/Normal/N_171'
sample=open('%s/%s_sample.txt'%(maindir,normal),'r')
sample_list=set([line.rstrip('\n') for line in sample.readlines()])
os.system("rm %s/result/%s_result.txt"%(maindir,normal))

for i in sample_list:
	os.system('cat %s/%s/star-fusion.fusion_predictions.abridged.tsv > %s/result/%s.txt'%(maindir,i,maindir,i))
	os.system('ls %s/result/%s.txt >> %s/result/%s_result.txt'%(maindir,i,maindir,normal))
	os.system('cat %s/result/%s.txt >> %s/result/%s_result.txt'%(maindir,i,maindir,normal))

os.system("sed '/#/d' %s/result/%s_result.txt > %s/result/1.txt"%(maindir,normal,maindir))
os.system("sed 's/\/storage2\/Project\/TCGA_fusion\/STARFusion\/Normal\/%s\/result\///g' %s/result/1.txt > %s/result/%s_result.txt"%(normal,maindir,maindir,normal))

result = open("%s/result/%s_result.txt"%(maindir,normal),'r') 
result_lines = result.readlines()
sample = re.compile(sample_pattern)
result_table = open("%s/result/%s_result_table.txt"%(maindir,normal),'w')
#sample_ID   fusion   JunctionReadCount   SpanningFragCount
for i in result_lines:
	sample_check = sample.match(i)
	if sample_check :
		sample_ID = i.rstrip(".txt\n")
	else :
		result_table.write("%s\t%s\t%s"%(normal,sample_ID,i))

result.close()
result_table.close()

result_T = open("%s/result/%s_result_table.txt"%(maindir,normal),'r') 
result_Tlines = result_T.readlines()
filtered_result = open("%s/result/%s_filtered_result.txt"%(maindir,normal),'w')
for i in result_Tlines:
	seed = int(i.split('\t')[3])
	spanning = int(i.split('\t')[4].rstrip('\n'))
	if (seed == 1 & spanning >= 2) | seed >= 2:
		filtered_result.write("%s"%i)

result_T.close()
filtered_result.close()


filtered_result = open("%s/result/%s_filtered_result.txt"%(maindir,normal),'r')
Flines = filtered_result.readlines()
para_filtered_result = open("%s/result/%s_filtered_paralog_result.txt"%(maindir,normal),'w')

for i in Flines:
	fusion = i.split('\t')[2]
	H_gene = fusion.split('--')[0]
	T_gene = fusion.split('--')[1]
	if len(H_gene) <= 3 | len(T_gene) <= 3:
		if H_gene[0:2] != T_gene[0:2]:
			para_filtered_result.write("%s"%i)
	elif H_gene[0:3] != T_gene[0:3]:
		para_filtered_result.write("%s"%i)

filtered_result.close()
para_filtered_result.close()

para_filtered_result = open("%s/result/%s_filtered_paralog_result.txt"%(maindir,normal),'r')
Plines = para_filtered_result.readlines()
para_IG_filtered_result = open("%s/result/%s_filtered_paralog_IG_result.txt"%(maindir,normal),'w')

for i in Plines:
	fusion = i.split('\t')[2]
	H_gene = fusion.split('--')[0]
	T_gene = fusion.split('--')[1]
	if H_gene.find("-") == -1 & T_gene.find("-") == -1 :
		if H_gene.find(".") == -1 & T_gene.find(".") == -1 :
			para_IG_filtered_result.write("%s"%i)

para_filtered_result.close()
para_IG_filtered_result.close()

os.system("cut -f3 %s/result/%s_filtered_paralog_IG_result.txt | sort | uniq -dc | sort > %s/result/%s_duplicate.txt"%(maindir,normal,maindir,normal))
