import numpy as np

def read_string(fname,delimiter=' ',comment='#'):
	f = open(fname,'r')
	first_time = True
	while True:
		line = f.readline()
		if line == '':
			break

		if line[0] != comment:
			elems = line.split()
			elems[-1] = elems[-1][:-1]
			if first_time:
				all_elems = np.array(elems)
				first_time = False
			else:
				all_elems = np.vstack((all_elems, np.array(elems)))
	return all_elems.transpose()

filename = raw_input('\t Filename of the file to format?')
sn = raw_input('\t Target signal-to-noise of the objects (default is 20, press enter to get the default)? ')
if sn == '':
    sn = '20'

name, ra, dec, mag, P, prio  = read_string(filename)

f = open('formatted_'+filename,'w')

f.write('type\tccdros\tsn\tnoprog\trefnocod\talphacat\tdeltacat\tequicat\tmv\n')
f.write('------\n')

for i in range(len(name)):
	f.write('COR_OBFP\tslow\t'+sn+'\t1\t'+name[i]+'\t'+\
			ra[i]+'\t'+dec[i]+'\t2000.0\t'+mag[i]+'\n')
f.close()
