import os,glob

main_dir = '/storage2/Project/TCGA_fusion/Raw_TCGA'
cancer = 'BRCA' 
jason = 'files.2018-06-04.json' 
tsv = 'BRCA.tsv'

os.system("cd %s/%s/fasta"%(main_dir,cancer))
os.system("ls -R | grep './\|.tar.gz' | sed '/logs/d'> %s_tar_list.txt"%cancer)

Fdict = {}
file = open("%s/%s/fasta/%s_tar_list.txt"%(main_dir,cancer,cancer))
for line in file:
    try:
        if './' in line:
            folder = line.split('./')[1].split(':')[0]
        elif '.tar.gz' in line:
            tar = line.split('\n')[0]
        Fdict[folder] = tar
    except:
        print line


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

Bdict={}
gdc = open("%s/%s/%s"%(main_dir,cancer,tsv))
gdc_lines = gdc.readlines()
for line in gdc_lines:
    try: 
        uuid = line.split('\t')[2]
        bar = line.split('\t')[3]
        Bdict[uuid]=bar
    except:
        print line

for folder in Fdict.keys():
    print "%s\t%s"%(folder, Bdict[Gdict[Fdict[folder]]])
####
####simbolic link
lndir="%s/%s/ln_fasta"%(main_dir,cancer)

for folder in Fdict.keys():
    lninf=glob.glob("%s/%s/fasta/%s/*.fastq"%(main_dir,cancer,folder))
    for lnline in lninf:
        rfile="_".join(lnline.split("/")[-1].split("_")[:-1])
        inputdir1="/".join(lnline.split("/")[:-1])
        inputfile1="%s/%s_1.fastq"%(inputdir1, rfile); inputfile2="%s/%s_2.fastq"%(inputdi
r1, rfile)
        EXfile1=os.path.isfile(inputfile1); EXfile2=os.path.isfile(inputfile2)
        print inputfile1,inputfile2
        if EXfile1 == False or EXfile2 ==False:
            print "no file",inputfile1,inputfile1
            continue
        os.system("ln -s %s %s/%s%s"%(inputfile1,lndir,Bdict[Gdict[Fdict[folder]]],"_1.fastq"))
        os.system("ln -s %s %s/%s%s"%(inputfile2,lndir,Bdict[Gdict[Fdict[folder]]],"_2.fastq"))
