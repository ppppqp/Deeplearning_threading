#!/mnt/home/wuyunqiq/bin/python/anaconda2/bin/python
import os
import sys
import commands
import getopt
import shutil
import random

################################################### set home directory ############################################
os.environ['Met_predictor_HOME']='/nfs/amino-home/panqp/protein_stru/repo/met_predictor/Met-Predictor_12302019';
HOME=os.environ.get('Met_predictor_HOME');
PtmGetFeatures=HOME+'/bin/PtmGetFeatures.exe';
PtmGetFeatures2=HOME+'/bin/addstructurefeature.py';
LibsvmDir=HOME+'/lib/libsvm-3.14';
mono=HOME+'/lib/mono/bin/mono';
mypython='/usr/bin/python'
################################################### check file and set variable ###################################
isstructure = '1';
isscale = '1';
keepfea = 0;
 
opts,args = getopt.getopt(sys.argv[1:],"h:help:i:o:t:s:r:d:");
for op, value in opts:
    if op == "-h" or op == "-help":
        print 'usage: [python] Run_Metpredictor.py -i input_fasta -o outfile -t type [-s structure -r isscale]\n';
	print '-i input_fasta : use absolute path\n'
        print '-o outfile : use absolute path\n'
        print '-t type=K or R\n';
	print '-s structure: 1 for adding structure; 0 for not\n'
	print '-r isscale: 1 for scale; 0 for not\n' 
        print '-d keep features file in this directory\n'
        print 'example: ./Run_Metpredictor.py -i /home/Met-predictor/example/P0CX53.fasta -o /home/Met-predictor/example/P0CX53.out -t R -s 1 -r 1 -d ./P0CX53_fea \n';
        exit(0);
    elif op == "-i":
        inputfasta = os.path.abspath(value);
    elif op == "-o":
        outfile = os.path.abspath(value);
    elif op == "-t":
        residue = value;
    elif op == "-s":
        isstructure = value;
    elif op == "-r":
        isscale = value;
    elif op == "-d":
        feadir = os.path.abspath(value);
        os.system("mkdir -p "+feadir)
        keepfea = 1;

(Status,fastabasename)=commands.getstatusoutput('basename '+inputfasta+' .fasta');

AAindex=HOME+'/data/AAindex.dat';
KNNsampleFile=HOME+'/data/sample/'+residue+'.sample';

if residue == 'K':
    types = ['K','KMONO','KDI','KTRI'];
if residue == 'R':
    types = ['R','RMONO','RDI'];


SelectionFile=[];
RangeFile=[];
ModelFile=[];
for node in types:
    temp = node + '.' + isscale + '.' + isstructure;
    SelectionFile.append(HOME+'/data/model/'+temp+'/'+node+'.list');
    RangeFile.append(HOME+'/data/model/'+temp+'/'+node+'.range');
    ModelFile.append(HOME+'/data/model/'+temp+'/'+node+'.model');

################################### run soft  ####################################################
(Status,CurrentDir)=commands.getstatusoutput('pwd');
#FeatureDir=HOME+'/features/';
FeatureDir='/tmp/'+os.getenv('USER')+'/'+fastabasename+'-'+str(random.randint(999999999,9999999999))
print FeatureDir
os.system("rm -f "+FeatureDir+"/*")
if not os.path.exists(FeatureDir):
    os.makedirs(FeatureDir);
os.chdir(FeatureDir);


print 'calculate features...\n';


os.system('cp '+HOME+'/scripts/GetFeature.sh ./GetFeature.sh')
(Status,Info)=commands.getstatusoutput('./GetFeature.sh '+inputfasta);
print Info+'\n';

if isstructure == '1':
    print 'calculate Structure features...\n';
    os.system('cp '+HOME+'/scripts/GetStructureFeature.sh ./GetStructureFeature.sh')
    (Status,Info)=commands.getstatusoutput('./GetStructureFeature.sh '+inputfasta);
    print Info+'\n';
    os.system("mv "+fastabasename+"_features "+fastabasename+"-features")

#os.system("cp -r "+feadir+"/* "+FeatureDir);

if keepfea==1:
    os.system("cp -r "+FeatureDir+"/* "+feadir);
os.chdir(CurrentDir);
################################### calculate features ####################################################
windows=[];
if residue == 'K' and isstructure == '1':
    windows = ['17','23','15','23'];
if residue == 'K' and isstructure == '0':
    windows = ['17','23','17','21'];
if residue == 'R' and isstructure == '1':
    windows = ['17','17','21'];
if residue == 'R' and isstructure == '0':
    windows = ['17','17','19'];

for ii in range(0,len(windows)):
    window = windows[ii];
    
    libsvmFileold='./'+types[ii]+'.libsvm.old';
    (Status,Info)=commands.getstatusoutput(mono+' '+PtmGetFeatures+' '+inputfasta+' '+AAindex+' '+FeatureDir+' '+KNNsampleFile+' '+window+' '+residue+' '+libsvmFileold);
    print Info+'\n';

    if isstructure == '1':
        libsvmFile='./'+types[ii]+'.libsvm.new';
        (Status,Info)=commands.getstatusoutput(mypython+' '+PtmGetFeatures2+' '+FeatureDir+'/'+fastabasename+'-features/'+' '+residue+' '+window+' '+fastabasename+' '+inputfasta+' '+libsvmFileold+' '+libsvmFile);
        print Info+'\n';
        shutil.copyfile(libsvmFile,'./'+types[ii]+'.libsvm');
    else:
        shutil.copyfile(libsvmFileold,'./'+types[ii]+'.libsvm');

#os.chdir(CurrentDir);

####################################### clean features directory ####################################################
#if keepfea==1:
#    os.system("cp -r "+FeatureDir+"/* "+feadir);
os.system('rm -rf '+FeatureDir);


######################################### select features #################################################
print 'select featurs...\n';

selector=HOME+'/bin/PTM.selectcor.py';
for ii in range(0,len(types)):
    libsvmFile='./'+types[ii]+'.libsvm';
    outlibsvmFile='./'+types[ii]+'.out.libsvm';
    (Status,Info)=commands.getstatusoutput(mypython+' '+selector+' '+libsvmFile+' '+SelectionFile[ii]+' '+outlibsvmFile);


######################################### scale features #################################################
print 'scale libsvm file\n';
svmscale=LibsvmDir+'/svm-scale';
for ii in range(0,len(types)):
    outlibsvmFile='./'+types[ii]+'.out.libsvm';
    scalelibsvmFile='./'+types[ii]+'.scale.libsvm';
    if isscale == '1':
        (Status,Info)=commands.getstatusoutput(svmscale+' -r '+RangeFile[ii]+' '+outlibsvmFile+' >'+scalelibsvmFile);
    else:
        shutil.copyfile(outlibsvmFile,scalelibsvmFile);


############################################# run libsvm predict ####################################################
print 'predicting...\n';
svmpredictor=LibsvmDir+'/svm-predict';
for ii in range(0,len(types)):
    scalelibsvmFile='./'+types[ii]+'.scale.libsvm';
    svmoutFile='./'+types[ii]+'.svm.out';
    (Status,Info)=commands.getstatusoutput(svmpredictor+' -b 1 '+scalelibsvmFile+' '+ModelFile[ii]+' '+svmoutFile);


########### get index  ###############
sequence = '';
ifr = open(inputfasta,'r');
ifr.readline();
while 1:
    line = ifr.readline().strip();
    if not line:
        break;
    sequence = sequence + line;
ifr.close();

indexlist=[];
for ii in range(0,len(sequence)):
    if sequence[ii] == residue:
        indexlist.append(ii+1);


############################################# print outfile ####################################################
linelabel = [];
num = 0;
ifr = open('./'+types[0]+'.svm.out','r');
ifr.readline();
while 1:
    line = ifr.readline().strip();
    if not line:
        break;
    if line.split()[0] == '1':
        linelabel.append(num);
    num = num + 1;
ifr.close();


outlines = [];
for ii in range(0,len(types)):
    temp = [];
    svmoutFile='./'+types[ii]+'.svm.out';
    ifr = open(svmoutFile,'r');
    ifr.readline();
    while 1:
        line = ifr.readline().strip();
        if not line:
            break;
        temp.append(line);
    ifr.close();
    outlines.append(temp);


ofw = open(outfile,'w');
ofw.write('SeqIndex Label Pos_prob Neg_prob MONOLabel MONOPos_prob MONONeg_prob DILabel DIPos_prob DINeg_prob TRILabel TRIPos_prob TRINeg_prob\n');
for jj in range(0,len(indexlist)):
    ofw.write(str(indexlist[jj])+' '+outlines[0][jj]);
    if jj in linelabel:
        for ii in range(1,len(types)):
            ofw.write(' '+outlines[ii][jj]);
    ofw.write('\n');
ofw.close();

################################################## clean files ########################################################
for ii in range(0,len(types)):
    libsvmFileold='./'+types[ii]+'.libsvm.old';
    libsvmFilenew='./'+types[ii]+'.libsvm.new';
    libsvmFile='./'+types[ii]+'.libsvm';
    outlibsvmFile='./'+types[ii]+'.out.libsvm';
    scalelibsvmFile='./'+types[ii]+'.scale.libsvm';
    svmoutFile='./'+types[ii]+'.svm.out';
    if isstructure == '1':
        os.system('rm '+libsvmFileold+' '+libsvmFilenew+' '+libsvmFile+' '+outlibsvmFile+' '+scalelibsvmFile+' '+svmoutFile);
    else:
        os.system('rm '+libsvmFileold+' '+libsvmFile+' '+outlibsvmFile+' '+scalelibsvmFile+' '+svmoutFile);

