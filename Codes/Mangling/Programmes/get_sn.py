## Import SN ###
data_sn=open('%s_sifto.txt'%sn_name,'r')
lines = data_sn.readlines()
fields = lines[0].split()
if len(fields) != 5: 
	print('first line of %s must have 6 fields:  name, redshift, zcmb, RA, DEC'%sn_name)
name = fields[0]
try:
	z_sn,z_sn_cmb,ra,decl= map(float, fields[1:5])
except:
	print("z, ra, and dec must be floats  (ra/dec in decimal degrees)")


lines = lines[1:]
this_filter = None
bands=[]
MJD = {}
mags = {}
emags = {}

for line in lines:
	if line[0] == "#":  continue
	if line.find('filter') >= 0:
		this_filter = line.split()[1]
		bands.append(this_filter)
		MJD[this_filter] = []
		mags[this_filter] = []
		emags[this_filter] = []
	elif this_filter is not None:
		try:
			t,m,em = map(float, str.split(str.strip(line)))
		except:
			print('Bad format in line:\n %s' % (line))
		MJD[this_filter].append(t)
		mags[this_filter].append(m)
		emags[this_filter].append(em)

for f in MJD:
	MJD[f] = array(MJD[f])
	mags[f] = array(mags[f])
	emags[f] = array(emags[f])
