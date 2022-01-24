#!usr/bin/python
import sys;
import os;
# import commands;
import math;

def get_status_output(*args, **kwargs):
    out = subprocess.run(*args, **kwargs)
    return out.stdout

    return p.returncode, stdout, stderr
def getfragment(fastaname,fragmentdir):
    begins = [];
    ends = [];
    outstr = get_status_output('ls '+'/nfs/amino-home/panqp/protein_stru/repo/met_predictor/Met-Predictor_12302019/example/example_features/P0CX53_fasta/P0CX53_s*.fasta', shell=True, capture_output=True);
    print(outstr)
    #print(fragmentdir)
    #print(outstr)
    for items in outstr.split('\n'):
        #print(items)
        begins.append(int(items.split('.')[0].split('_')[1][1:]));
        ends.append(int(items.split('.')[0].split('_')[2][1:]));
    return begins,ends;

def getchops(chopsfile,site,window):
    chops = [];
    # begin = int(chopsfile.split('/')[-1].split('.')[0].split('_')[1][1:]);
    begin = 1
    halfwindow = int((window - 1)/ 2)
    newsite = site - begin + 1 + halfwindow;

    outlist = [];
    ifr = open(chopsfile,'r');
    for line in ifr:
        if line.strip()[13:15] == 'CA':
            x = float(line.strip()[60:66].strip())/100;
            x = min(1.0,x);
            outlist.append(str(round(x,3))); 
    ifr.close();
    outlist = ['0' for x in range(0,halfwindow)] + outlist + ['0' for x in range(0,halfwindow)];

    chops = outlist[newsite-halfwindow-1:newsite+halfwindow];
    return chops;

def gethse(hsefile,site,window):
    hse = [];
    # begin = int(hsefile.split('/')[-1].split('.')[0].split('_')[1][1:]);
    begin = 1
    halfwindow = int((window - 1)/2);
    newsite = site - begin + 1 + halfwindow;

    outlist = [];
    ifr = open(hsefile,'r');
    for line in ifr:
        if line.strip()[13:15] == 'CA':
            x = float(line.strip()[60:66].strip())/100;
            x = min(1.0,x);
            outlist.append(str(round(x,3))); 
    ifr.close();
    outlist = ['0' for x in range(0,halfwindow)] + outlist + ['0' for x in range(0,halfwindow)];

    hse = outlist[newsite-halfwindow-1:newsite+halfwindow];
    return hse;

def getkth(kthfile,site,KR):
    kth = '0';
    
    # begin = int(kthfile.split('/')[-1].split('.')[0].split('_')[1][1:]);
    begin = 1
    newsite = site - begin + 1;
    
    atoms = [];
    if KR == 'K':
        atoms = ['NZ '];
    if KR == 'R':
        atoms = ['NH1','NH2'];

    tmp = []
    ifr = open(kthfile,'r');
    for line in ifr:
        if line.strip()[13:16] in atoms and line.strip()[22:26].strip() == str(newsite):
            x = float(line.strip()[60:66].strip())/100;
            x = min(1.0,x);
            tmp.append(x);
    ifr.close();
    
    kth = str(sum(tmp)/float(len(tmp)));
    return kth;

def getl1depth(l1file,site,KR):
    l1 = [];
        
    # begin = int(l1file.split('/')[-1].split('.')[0].split('_')[1][1:]);
    begin = 1
    newsite = site - begin + 1;
    
    atoms = [];
    if KR == 'K':
        atoms = ['NZ '];
    if KR == 'R':
        atoms = ['NH1','NH2'];

    tmp = []
    ifr = open(l1file,'r');
    for line in ifr:
        if line.strip()[13:16] in atoms and line.strip()[22:26].strip() == str(newsite):
            x = float(line.strip()[54:60].strip());
            tmp.append(round(x,3));
    ifr.close();

    l1.append(str(sum(tmp)/float(len(tmp))));
    return l1;

def getresiduedepth(residuefile,site,KR):
    residuedepth = '0';
    # 
    # begin = int(residuefile.split('/')[-1].split('.')[0].split('_')[1][1:]);
    begin = 1
    newsite = site - begin + 1;
    
    atoms = [];
    if KR == 'K':
        atoms = ['NZ '];
    if KR == 'R':
        atoms = ['NH1','NH2'];

    tmp = []
    ifr = open(residuefile,'r');
    for line in ifr:
        if line.strip()[13:16] in atoms and line.strip()[22:26].strip() == str(newsite):
            x = float(line.strip()[60:66].strip())/20;
            x = min(1.0,x);
            tmp.append(x);
    ifr.close();
    
    residuedepth = str(sum(tmp)/float(len(tmp)));
    return residuedepth;

def getrsa(rsafile,site,window):
    rsa = [];
    
    # begin = int(rsafile.split('/')[-1].split('.')[0].split('_')[1][1:]);
    begin = 1
    halfwindow = int((window - 1)/2);
    newsite = site - begin + 1 + halfwindow;

    outlist2 = [];
    ifr = open(rsafile,'r');
    for line in ifr:
        if line.startswith('RES'):
            x = float(line.strip()[35:41].strip())/100;
            x = min(1.0,x);
            outlist2.append(str(x));
    ifr.close();
    outlist2 = ['0' for x in range(0,halfwindow)] + outlist2 + ['0' for x in range(0,halfwindow)];

    rsa = outlist2[newsite-halfwindow-1:newsite+halfwindow];
    return rsa;

def getglobalkth(globalkthfile,site,KR):
    globalkth = '0';
    
    atoms = [];
    if KR == 'K':
        atoms = ['NZ '];
    if KR == 'R':
        atoms = ['NH1','NH2'];

    
    if os.path.exists(globalkthfile):
        tmp = [];
        ifr = open(globalkthfile,'r');
        for line in ifr:
            if line.strip()[13:16] in atoms and line.strip()[22:26].strip() == str(site):
                x = float(line.strip()[60:66].strip())/100;
                x = min(1.0,x);
                tmp.append(x);
        ifr.close();
        globalkth = str(sum(tmp)/float(len(tmp)));
    else:
        print (globalkthfile)
    
    
    return globalkth;

def gettmscore(tmscorefile):
    tmscore = '0';
    ifr = open(tmscorefile,'r');
    for line in ifr:
        if line.startswith('REMARK'):
            tmscore = line.strip().split()[2];
            break;
    ifr.close();
    return tmscore;
    
def writefile(fastaname,KR,fragmentdir,chopsdir,hsedir,kthdir,l1dir,residuedir,rsadir,globalkthdir,tmscoredir,sitelist,window,index):
    outfilelist = [];
    # (begins,ends) = getfragment(fastaname,fragmentdir);
    for node in range(0,len(sitelist)):
        site = sitelist[node];
        ########################
        # begin = '';
        # end = '';
        # for ii in range(0,len(begins)):
        #     if site >= begins[ii] and site <= ends[ii]:
        #         begin = str(begins[ii]);
        #         end = str(ends[ii]);
        ########################
        chopsfile = chopsdir + fastaname + '.sch.pdb';
        chops = getchops(chopsfile,site,window);

        hselist = [];
        for ii in range(5,31):
            for hsetype in ['_HSEAD','_HSEAU','_HSEBD','_HSEBU']:
                hsefile = hsedir + str(ii) + '/' + fastaname + hsetype + str(ii) + '.pdb';
                hse = gethse(hsefile,site,window);
                hselist = hselist + hse;

        kthfile = kthdir + fastaname  + '.ch.pdb';
        kth = getkth(kthfile,site,KR);

        l1list = [];
        for ii in range(5,31):
            l1file = l1dir + str(ii) + '/' + fastaname + '_atom-atom_local' + str(ii) + '.pdb';
            l1 = getl1depth(l1file,site,KR);
            l1list = l1list + l1;

        residuefile = residuedir + fastaname + '-atomic_depth.pdb';
        residuedepth = getresiduedepth(residuefile,site,KR);

        rsafile = rsadir + fastaname +'.rsa';
        rsa = getrsa(rsafile,site,window);

        globalkthfile = globalkthdir + fastaname.replace('-','_') +'.ch.pdb';
        globalkth = getglobalkth(globalkthfile,site,KR);

        tmscorefile = tmscoredir + fastaname + '.pdb';
        tmscore = gettmscore(tmscorefile);
        
        ########################
        outfileline = '';
        for ii in chops:
            outfileline = outfileline + ' ' + str(index) + ':' + ii;
            index = index  + 1;
        for ii in hselist:
            outfileline = outfileline + ' ' + str(index) + ':' + ii;
            index = index  + 1;
        outfileline = outfileline + ' ' + str(index) + ':' + kth;
        index = index  + 1;
        for ii in l1list:
            outfileline = outfileline + ' ' + str(index) + ':' + ii;
            index = index  + 1;
        outfileline = outfileline + ' ' + str(index) + ':' + residuedepth;
        index = index  + 1;
        for ii in rsa:
            outfileline = outfileline + ' ' + str(index) + ':' + ii;
            index = index  + 1;
        outfileline = outfileline + ' ' + str(index) + ':' + globalkth;
        index = index  + 1;
        outfileline = outfileline + ' ' + str(index) + ':' + tmscore;
        index = index  + 1;
        outfilelist.append(outfileline);
        
    return outfilelist;
    
if __name__ == '__main__':    
    #fragmentdir = '/home/wuyun/Met/dataset/fragmentfasta/';
    #chopsdir = '/home/wuyun/Met/dataset/fragmentchops/';
    #hsedir = '/home/wuyun/Met/dataset/fragmenthse/';
    #kthdir = '/home/wuyun/Met/dataset/fragmentkthCH/';
    #l1dir = '/home/wuyun/Met/dataset/fragmentL1depth/';
    #residuedir = '/home/wuyun/Met/dataset/fragmentresiduedepth/';
    #rsadir = '/home/wuyun/Met/dataset/fragmentrsa/';
    #globalkthdir = '/home/wuyun/Met/dataset/pdbkthCH/';
    #fastaname = 'A1X283';#P20810-5
    #sitefasta = '/home/wuyun/Met/seqtrain/site_format/Bench_K/'+fastaname+'.fasta';
    #tmscoredir = '';
    #sitelist = [619];
    #window = 9;
    #KR = 'K';
    #writefile(fastaname,KR,fragmentdir,chopsdir,hsedir,kthdir,l1dir,residuedir,rsadir,globalkthdir,tmscoredir,sitelist,window,1);



    featuredir = sys.argv[1]+'/';
    KR = sys.argv[2];
    window = int(sys.argv[3]);
    fastaname = sys.argv[4];
    fastafile = sys.argv[5];
    seqlibsvmfile = sys.argv[6];
    outlibsvmfile = sys.argv[7];
    

    
    sequence = '';
    ifr = open(fastafile,'r');
    ifr.readline();
    while 1:
        line = ifr.readline().strip();
        if not line:
            break;
        sequence = sequence + line;
    ifr.close();
    sitelist = [];
    for ii in range(0,len(sequence)):
        if sequence[ii] == KR:
            sitelist.append(ii+1);


    num = 0;
    ofw = open(outlibsvmfile,'w');
    sfr = open(seqlibsvmfile,'r');
    for line in sfr:
        site = [sitelist[num]];
        index = int(line.strip().split()[-1].split(':')[0])+1;
        outfilelist = writefile(fastaname,KR,featuredir,featuredir,featuredir,featuredir,featuredir,featuredir,featuredir,featuredir,featuredir,site,window,index);
        ofw.write(line.strip()+' '+outfilelist[0].strip()+'\n');
        num = num+1;
    sfr.close();
    ofw.close();
                

    
    
