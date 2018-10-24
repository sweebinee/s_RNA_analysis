

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



sample_list=['00c8bbd7-dd91-43bd-ab18-5dc37911feb0','03321e5d-700e-4877-9857-40d7310a3936','05d92117-aaf6-4503-9b55-822af0edfe06','0675d290-0bd7-4114-9fd5-924cfa82fc27','078240cf-327e-4c53-bdc8-a3dfaa55ebe4','07931666-ca61-4446-99ce-6cd9f9b71133','085c85bb-6c3c-4c60-9d35-35d8b57b894f','08fc4125-35e7-43b9-ac3c-2883bf879240','0900a4c2-8e3d-48be-978b-5507dbcad22c','0bf47a03-0532-4f5f-bde2-4b801bd18ccb','0df186bf-6a8c-4777-835a-80b621e46ad6','0f1e80db-5d42-4ec1-82ba-24bc4361e05d','0fd07ec8-3021-4158-b7a2-bd90d2c56ae3','10547eba-1f81-4861-8e9f-3f01a290a44f','11bc2c5d-1ceb-4542-a527-4ff843e3c04d','126d9772-ad58-48c4-b5b8-def3bec1ea03','12ab5ebe-d9e0-41e8-b555-291aa6bd9e82','12c29388-dadf-4291-8226-465641fd11ec','1395d9ec-8d5f-4b21-a709-a36b827f27db','1477fb6b-9090-48ca-b3a8-436f893e80e2','15730dab-048f-4606-b5f9-1c013a79f426','16f0aa6b-3e61-4173-b8f6-994911070d17','18780597-255a-42f2-b06b-f986d7bcbcf9','18da7680-c797-46bb-a926-e99fc5a67abb','19fd6073-3ef3-4ab2-bc00-bb6873035f5f','1ba41b90-9257-4fac-99a1-083d9ab27024','1be7635d-6bad-4368-ae22-9c4df63a77ac','1cd1feec-67fe-4c00-9ab8-4dc5e01375ac','1dd9a7cf-90ef-4a1b-a3ca-09ddbf9661e3','1f669df8-b53b-4b7a-a170-135d7603c555']
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

sample_list=['5fc0927e-070b-463d-ad14-892d2beba700']
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





    