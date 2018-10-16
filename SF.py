

import os,glob,datetime
cancer="STAD"
Gdict={}
inf=open("/storage2/Project/TCGA_fusion/Raw_TCGA/%s/STADfasta.xml"%cancer)
for line in inf:
    try:
        if "<analysis_id>" in line:
            uuid = line.split("<analysis_id>")[1].split("</analysis_id>")[0]
        elif "<legacy_sample_id>" in line:
            bar = line.split("<legacy_sample_id>")[1].split("</legacy_sample_id>")[0]
        Gdict[uuid]=bar
    except:
        print line


lninf=glob.glob("/storage2/Project/TCGA_fusion/Raw_TCGA/BRCA/fasta/00a4e56f-c9f4-4c3d-b795-9b468de83223/*.fastq")

lndir="/storage2/Project/TCGA_fusion/Raw_TCGA/%s/ln_fasta"%(cancer)

runlist=[]
for lnline in lninf:
    rfile="_".join(lnline.split("/")[-1].split("_")[:-1])
    inputdir1="/".join(lnline.split("/")[:-1])
    inputfile1="%s/%s_1.fastq"%(inputdir1, rfile); inputfile2="%s/%s_2.fastq"%(inputdir1, rfile)
    EXfile1=os.path.isfile(inputfile1); EXfile2=os.path.isfile(inputfile2)
    print inputfile1,inputfile2
    if EXfile1 == False or EXfile2 ==False:
        print "no file",inputfile1,inputfile1
        continue
    uuid=lnline.split("/")[-2]
    os.system("ln -s %s %s/%s%s"%(inputfile1,lndir,Gdict[uuid],"_1.fastq"))
    os.system("ln -s %s %s/%s%s"%(inputfile2,lndir,Gdict[uuid],"_2.fastq"))
    if Gdict[uuid] not in runlist:
        runlist.append(Gdict[uuid])



sample_list=['0675d290-0bd7-4114-9fd5-924cfa82fc27']
for sample in sample_list:
    lninf=glob.glob("/storage2/Project/TCGA_fusion/Raw_TCGA/%s/fasta/%s/*.fastq"%(cancer,sample))
    lndir="/storage2/Project/TCGA_fusion/Raw_TCGA/%s/ln_fasta"%(cancer)
    rfile="_".join(lninf[0].split("/")[-1].split("_")[:-1])
    inputdir1="/".join(lninf[0].split("/")[:-1])
    inputfile1="%s/%s_1.fastq"%(inputdir1, rfile); inputfile2="%s/%s_2.fastq"%(inputdir1, rfile)
    EXfile1=os.path.isfile(inputfile1); EXfile2=os.path.isfile(inputfile2)
    if EXfile1 == False or EXfile2 ==False:
        print "no file",inputfile1,inputfile1
    uuid=lninf[0].split("/")[-2]
    os.system("ln -s %s %s/%s%s"%(inputfile1,lndir,Gdict[uuid],"_1.fastq"))
    os.system("ln -s %s %s/%s%s"%(inputfile2,lndir,Gdict[uuid],"_2.fastq"))
    print Gdict[uuid]



    