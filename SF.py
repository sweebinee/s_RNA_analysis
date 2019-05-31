

import os,glob,datetime
cancer="STAD"
Gdict={}
inf=open("/storage2/Project/TCGA_fusion/Raw_TCGA/%s/STADbam_20141029.xml"%cancer)
for line in inf:
    try:
        if "<analysis_id>" in line:
            uuid = line.split("<analysis_id>")[1].split("</analysis_id>")[0]
        elif "<legacy_sample_id>" in line:
            bar = line.split("<legacy_sample_id>")[1].split("</legacy_sample_id>")[0]
        Gdict[uuid]=bar
    except:
        print line


lninf=glob.glob("/storage2/Project/TCGA_fusion/Raw_TCGA/BRCA/fasta/05d92117-aaf6-4503-9b55-822af0edfe06/*.fastq")

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



sample_list=['02fd3ddb-623b-41af-8edf-9cb95c18bb89']
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


#######################################################################################
#STAD fasta gz file has sample_id

sample_list=['12ab5ebe-d9e0-41e8-b555-291aa6bd9e82']
for sample in sample_list:
    lninf=glob.glob("/storage2/Project/TCGA_fusion/Raw_TCGA/%s/fasta/%s/*.tar"%(cancer,sample))
    #['/storage2/Project/TCGA_fusion/Raw_TCGA/STAD/fasta/00a4e56f-c9f4-4c3d-b795-9b468de83223/TCGA-D7-6818-01A-11R-1884-13_rnaseq_fastq.tar']
    lndir="/storage2/Project/TCGA_fusion/Raw_TCGA/%s/ln_fasta"%(cancer)
    rfile=lninf[0].split("/")[-1].split("_")[0]
    inputdir1="/".join(lninf[0].split("/")[:-1])
    inputfile1=glob.glob("%s/*_1.fastq"%inputdir1); inputfile2=glob.glob("%s/*_2.fastq"%inputdir1)
    EXfile1=os.path.isfile(inputfile1[0]); EXfile2=os.path.isfile(inputfile2[0])
    if EXfile1 == False or EXfile2 ==False:
        print "no file",inputfile1,inputfile1
    os.system("ln -s %s %s/%s%s"%(inputfile1[0],lndir,rfile,"_1.fastq"))
    os.system("ln -s %s %s/%s%s"%(inputfile2[0],lndir,rfile,"_2.fastq"))
    print rfile





    