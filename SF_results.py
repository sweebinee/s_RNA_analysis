#result RAW file collecting code (ver.STAR-Fusion)
ls -d TCGA* > PAAD_sample.txt

import os,re

cancer_list = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC']
cancer_list = ['KIRP','LGG','LIHC','LUAD','LUSC','MESO','OV','PCPG','PRAD','READ','SARC','SKCM']
cancer_list = ['STAD','THCA','THYM','UCEC','UCS','UVM','TGCT','PAAD']
#LAML 따로 

for cancer in cancer_list:
	#cancer = 'ACC'
	maindir = '/storage2/Project/TCGA_fusion/STARFusion/%s'%cancer
	sample = open("%s/%s_sample.txt"%(maindir,cancer),'r')
	sample_list = set([line.rstrip("\n") for line in sample.readlines()])
	os.system("rm %s/result/%s_result.txt"%(maindir,cancer))
	#i = 'TCGA-OR-A5LR-01A-11R-A29S-07'
	for i in sample_list:
		os.system("cp %s/%s/star-fusion.fusion_candidates.final %s/%s/star-fusion.fusion_predictions.abridged.tsv"%(maindir,i,maindir,i))
		os.system("cat %s/%s/star-fusion.fusion_predictions.abridged.tsv > %s/result/%s.txt"%(maindir,i,maindir,i))
		os.system("ls %s/result/%s.txt >> %s/result/%s_result.txt"%(maindir,i,maindir,cancer))
		os.system("cat %s/result/%s.txt >> %s/result/%s_result.txt"%(maindir,i,maindir,cancer))
	#
	os.system("sed '/#/d' %s/result/%s_result.txt > %s/result/1.txt"%(maindir,cancer,maindir))
	os.system("sed 's/\/storage2\/Project\/TCGA_fusion\/STARFusion\/%s\/result\///g' %s/result/1.txt > %s/result/%s_result.txt"%(cancer,maindir,maindir,cancer))
	#
	result = open("%s/result/%s_result.txt"%(maindir,cancer),'r') 
	result_lines = result.readlines()
	sample = re.compile("^TCGA")
	result_table = open("%s/result/%s_result_table.txt"%(maindir,cancer),'w')
	#sample_ID   fusion   JunctionReadCount   SpanningFragCount
	for i in result_lines:
		sample_check = sample.match(i)
		if sample_check :
			sample_ID = i.rstrip(".txt\n")
		else :
			result_table.write("%s\t%s\t%s"%(cancer,sample_ID,i))
	#
	result.close()
	result_table.close()
	#
	result_T = open("%s/result/%s_result_table.txt"%(maindir,cancer),'r') 
	result_Tlines = result_T.readlines()
	filtered_result = open("%s/result/%s_filtered_result.txt"%(maindir,cancer),'w')
	for i in result_Tlines:
		seed = int(i.split('\t')[3])
		spanning = int(i.split('\t')[4].rstrip('\n'))
		if (seed == 1 & spanning >= 2) | seed >= 2:
			filtered_result.write("%s"%i)
	#
	result_T.close()
	filtered_result.close()
	#
	os.system("cut -f2 %s/result/%s_filtered_result.txt | sort | uniq -dc | sort > %s/result/%s_duplicate.txt"%(maindir,cancer,maindir,cancer))
	#
	filtered_result = open("%s/result/%s_filtered_result.txt"%(maindir,cancer),'r')
	Flines = filtered_result.readlines()
	para_filtered_result = open("%s/result/%s_filtered_paralog_result.txt"%(maindir,cancer),'w')
	#
	for i in Flines:
		fusion = i.split('\t')[2]
		H_gene = fusion.split('--')[0]
		T_gene = fusion.split('--')[1]
		if len(H_gene) <= 3 | len(T_gene) <= 3:
			if H_gene[0:2] != T_gene[0:2]:
				para_filtered_result.write("%s"%i)
		elif H_gene[0:3] != T_gene[0:3]:
			para_filtered_result.write("%s"%i)
	#
	filtered_result.close()
	para_filtered_result.close()
	#
	para_filtered_result = open("%s/result/%s_filtered_paralog_result.txt"%(maindir,cancer),'r')
	Plines = para_filtered_result.readlines()
	para_IG_filtered_result = open("%s/result/%s_filtered_paralog_IG_result.txt"%(maindir,cancer),'w')
	#
	for i in Plines:
		fusion = i.split('\t')[2]
		H_gene = fusion.split('--')[0]
		T_gene = fusion.split('--')[1]
		if H_gene.find("-") == -1 & T_gene.find("-") == -1 :
			if H_gene.find(".") == -1 & T_gene.find(".") == -1 :
				para_IG_filtered_result.write("%s"%i)
	#
	para_filtered_result.close()
	para_IG_filtered_result.close()


#LAML
cancer = 'LAML'
maindir = '/storage2/Project/TCGA_fusion/STARFusion/%s'%cancer
sample = open("%s/%s_sample.txt"%(maindir,cancer),'r')
sample_list = set([line.rstrip("\n") for line in sample.readlines()])
os.system("rm %s/result/%s_result.txt"%(maindir,cancer))
#i = 'TCGA-OR-A5LR-01A-11R-A29S-07'
for i in sample_list:
	os.system("ls %s/result/%s.txt >> %s/result/%s_result.txt"%(maindir,i,maindir,cancer))
	os.system("cat %s/result/%s.txt >> %s/result/%s_result.txt"%(maindir,i,maindir,cancer))

os.system("sed '/#/d' %s/result/%s_result.txt > %s/result/1.txt"%(maindir,cancer,maindir))
os.system("sed 's/\/storage2\/Project\/TCGA_fusion\/STARFusion\/%s\/result\///g' %s/result/1.txt > %s/result/%s_result.txt"%(cancer,maindir,maindir,cancer))

result = open("%s/result/%s_result.txt"%(maindir,cancer),'r') 
result_lines = result.readlines()
sample = re.compile("^TCGA")
result_table = open("%s/result/%s_result_table.txt"%(maindir,cancer),'w')
#sample_ID   fusion   JunctionReadCount   SpanningFragCount
for i in result_lines:
	sample_check = sample.match(i)
	if sample_check :
		sample_ID = i.rstrip(".txt\n")
	else :
		result_table.write("%s\t%s\t%s"%(cancer,sample_ID,i))

result.close()
result_table.close()

result_T = open("%s/result/%s_result_table.txt"%(maindir,cancer),'r') 
result_Tlines = result_T.readlines()
filtered_result = open("%s/result/%s_filtered_result.txt"%(maindir,cancer),'w')
for i in result_Tlines:
	seed = int(i.split('\t')[3])
	spanning = int(i.split('\t')[4].rstrip('\n'))
	if (seed == 1 & spanning >= 2) | seed >= 2:
		filtered_result.write("%s"%i)

result_T.close()
filtered_result.close()

os.system("cut -f2 %s/result/%s_filtered_result.txt | sort | uniq -dc | sort > %s/result/%s_duplicate.txt"%(maindir,cancer,maindir,cancer))

filtered_result = open("%s/result/%s_filtered_result.txt"%(maindir,cancer),'r')
Flines = filtered_result.readlines()
para_filtered_result = open("%s/result/%s_filtered_paralog_result.txt"%(maindir,cancer),'w')

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

para_filtered_result = open("%s/result/%s_filtered_paralog_result.txt"%(maindir,cancer),'r')
Plines = para_filtered_result.readlines()
para_IG_filtered_result = open("%s/result/%s_filtered_paralog_IG_result.txt"%(maindir,cancer),'w')

for i in Plines:
	fusion = i.split('\t')[2]
	H_gene = fusion.split('--')[0]
	T_gene = fusion.split('--')[1]
	if H_gene.find("-") == -1 & T_gene.find("-") == -1 :
		if H_gene.find(".") == -1 & T_gene.find(".") == -1 :
			para_IG_filtered_result.write("%s"%i)

para_filtered_result.close()
para_IG_filtered_result.close()


_filtered_paralog_IG_result.txt

%s/result/%s_filtered_paralog_IG_result.txt

cp  /storage2/Project/TCGA_fusion/STARFusion/ACC/result/ACC_filtered_paralog_IG_result.txt  /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/BLCA/result/BLCA_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/BRCA/result/BRCA_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/CESC/result/CESC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/CHOL/result/CHOL_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/COAD/result/COAD_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/DLBC/result/DLBC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/ESCA/result/ESCA_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/GBM/result/GBM_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/HNSC/result/HNSC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/KICH/result/KICH_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/KIRC/result/KIRC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/KIRP/result/KIRP_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/LGG/result/LGG_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/LIHC/result/LIHC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/LUAD/result/LUAD_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/LUSC/result/LUSC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/MESO/result/MESO_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/OV/result/OV_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/PCPG/result/PCPG_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/PRAD/result/PRAD_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/READ/result/READ_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/SARC/result/SARC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/SKCM/result/SKCM_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/STAD/result/STAD_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/THCA/result/THCA_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/THYM/result/THYM_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/UCEC/result/UCEC_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/UCS/result/UCS_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/UVM/result/UVM_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/TGCT/result/TGCT_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/PAAD/result/PAAD_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
cp /storage2/Project/TCGA_fusion/STARFusion/LAML/result/LAML_filtered_paralog_IG_result.txt /storage2/Project/TCGA_fusion/STARFusion/result
############
#Tumor - Normal

cancer_list = set(['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PCPG','PRAD','READ','SARC','SKCM','STAD','THCA','THYM','UCEC','UCS','UVM','TGCT','PAAD'])

maindir = "/storage2/Project/TCGA_fusion/STARFusion/result"
		
#일단은 tumor랑 normal 나눠야함
for cancer in cancer_list:
	filtered_result = open('%s/%s_filtered_paralog_IG_result.txt'%(maindir,cancer),'r')
	Flines = filtered_result.readlines()
	Tumor_result = open('%s/%s_Tumor.txt'%(maindir,cancer),'w')
	Normal_result = open('%s/%s_Normal.txt'%(maindir,cancer),'w')
	for i in Flines:
		sample = i.split('\t')[1]
		if sample.find('-11A') != -1 or sample.find('-11B') != -1 :
			Normal_result.write('%s'%i)
		else:
			Tumor_result.write('%s'%i)
	filtered_result.close()
	Tumor_result.close()
	Normal_result.close()
	#normal에서 나온 gene-pair set()으로 정리
	normal_gene_pair = set([line.rstrip('\n').split('\t')[2] for line in open('%s/%s_Normal.txt'%(maindir,cancer),'r')])
	print(normal_gene_pair)
	#Tumor에서 나오면 제외함
	Tumor = open('%s/%s_Tumor.txt'%(maindir,cancer),'r')
	Tumor_line = Tumor.readlines()
	final = open('%s/%s_filtered_paralog_IG_TN_result.txt'%(maindir,cancer),'w')
	for i in Tumor_line :
		fusion = i.split('\t')[2]
		if fusion not in normal_gene_pair :
			final.write('%s'%i)
	Tumor.close()
	final.close()

#

#######
#(ver.FusionScan)
import os 

cancer = 'SKCM'ㅌ
maindir = '/storage/home/subin95/FusionScan_TCGA/%s/result'%cancer
sample = open("%s/%s_sample.txt"%(maindir,cancer),'r')
sample_list = set([line.rstrip("\n") for line in sample.readlines()])

for i in sample_list:
	os.system("cat %s/%s/%s_FusionScan_candidate_result.txt > %s/result/%s.txt"%(maindir,i,i,maindir,i))




import os 

cancer = 'KIRC'
main_dir = '/storage/home/subin/Fusionscan_TCGA/%s/result'%cancer

os.system("cd %s"%main_dir)
os.system("ls -d TCGA* > %s/%s.txt"%(main_dir,cancer))

sample = open("%s/%s.txt"%(main_dir,cancer))
sample_list = set([line.rstrip("\n") for line in sample.readlines()])

for i in sample_list:
	os.system("rm -r %s/%s/*_1"%(main_dir,i))
	os.system("rm -r %s/%s/*_2"%(main_dir,i))
	os.system("rm -r %s/%s/*_3"%(main_dir,i))
	os.system("rm -r %s/%s/*_4"%(main_dir,i))
	os.system("rm -r %s/%s/*_5"%(main_dir,i))
	os.system("rm -r %s/%s/*_6"%(main_dir,i))
	os.system("rm -r %s/%s/*_7"%(main_dir,i))
	os.system("rm -r %s/%s/*_8"%(main_dir,i))
	os.system("rm -r %s/%s/fa_dir"%(main_dir,i))
	os.system("rm -r %s/%s/genefa_dir"%(main_dir,i))
	os.system("rm -r %s/%s/log_dir"%(main_dir,i))
	os.system("rm -r %s/%s/pslx_dir"%(main_dir,i))
	os.system("rm -r %s/%s/rescue_dir"%(main_dir,i))
	os.system("rm -r %s/%s/temp_dir"%(main_dir,i))


################################################################3
#KOBIC result re_naming
#######################
$ ls > PAAD_file.txt

import os,glob

main_dir = '/storage2/Project/TCGA_fusion/STARFusion'
cancer = 'PAAD' #암종명
jason = 'files.2018-08-16.json' 
tsv = '%s.tsv'%cancer
manifest = open('%s/%s/gdc_manifest_20180816_120356.txt'%(main_dir,cancer),'r')
manifest.readline()
manifest_lines = manifest.readlines()

Mdict = {}
for i in manifest_lines:
	folder = i.split('\t')[0]
	file = i.split('\t')[1]
	Mdict[folder] = file

manifest.close()

Gdict = {}
gdc_legacy = open("%s/%s/%s"%(main_dir,cancer,jason))
for line in gdc_legacy:
    try:
        if 'file_name' in line:
            tar = line.split('"file_name": "')[1].split('",')[0]
        elif 'case_id' in line:
            uuid = line.split('"case_id": "')[1].split('"')[0]
        Gdict[tar]=uuid
    except:
        print line

gdc_legacy.close()

Bdict={}
gdc = open("%s/%s/%s"%(main_dir,cancer,tsv))
gdc_lines = gdc.readlines()
for line in gdc_lines:
    try: 
        uuid = line.split('\t')[1]
        bar = line.split('\t')[2]
        Bdict[uuid]=bar
    except:
        print line

gdc.close()

result = open('%s/%s/%s_file.txt'%(main_dir,cancer,cancer),'r')
result_list = set([line.rstrip("\n") for line in result.readlines()])
for i in result_list:
	os.system("mv %s/%s/%s %s/%s/%s"%(main_dir, cancer, i ,main_dir, cancer, Bdict[Gdict[Mdict[i]]]))













