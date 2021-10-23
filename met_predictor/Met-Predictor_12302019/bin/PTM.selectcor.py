#!/mnt/home/wuyunqiq/bin/python/anaconda2/bin/python
import sys;
import os;
import commands;

def selectcorfeature(featurefile,corfile,selectfeaturefile):
	corlist = [];
	cr = open(corfile,"r");
	for line in cr:
		if len(line) > 1:
			corlist.append(int(line.strip().split(' ')[0]));
	cr.close();

	sw = open(selectfeaturefile,"w");
	fr = open(featurefile,"r");
	for line in fr:
		if len(line) > 1:
			outstr = line.strip().split(' ')[0];
			for ii in range(1,len(corlist)+1):
				outstr = outstr+" "+str(ii)+":"+str(line.strip().split(' ')[corlist[ii-1]].split(':')[1]);
			sw.write(outstr+'\n');
	fr.close();
	sw.close();


if __name__ == '__main__':
	featurefile = sys.argv[1];
	corfile = sys.argv[2];
	selectfeaturefile = sys.argv[3];
	
	selectcorfeature(featurefile,corfile,selectfeaturefile);




