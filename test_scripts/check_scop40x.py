#!/usr/bin/python3

import re
import sys

sfonly = False
if len(sys.argv) > 1:
	assert len(sys.argv) == 2
	assert sys.argv[1] == "sf"
	sfonly = True

outdir = "../big_scop40x/"

platform = "linux"

ref_tup2time={(platform, 'fam', 'kappa'): 118, (platform, 'fam', 'all'): 203, (platform, 'sf', 'kappa'): 82, (platform, 'sf', 'all'): 148, (platform, 'fold', 'kappa'): 97, (platform, 'fold', 'all'): 161, ('linux', 'fam', 'kappa'): 95, ('linux', 'fam', 'all'): 161, ('linux', 'sf', 'kappa'): 96, ('linux', 'sf', 'all'): 170, ('linux', 'fold', 'kappa'): 96, ('linux', 'fold', 'all'): 166}
ref_tup2mem={(platform, 'fam', 'kappa'): '10.9', (platform, 'fam', 'all'): '7.91', (platform, 'sf', 'kappa'): '12.1', (platform, 'sf', 'all'): '9.19', (platform, 'fold', 'kappa'): '11.4', (platform, 'fold', 'all'): '8.40', ('linux', 'fam', 'kappa'): '8.41', ('linux', 'fam', 'all'): '6.49', ('linux', 'sf', 'kappa'): '9.67', ('linux', 'sf', 'all'): '7.76', ('linux', 'fold', 'kappa'): '8.88', ('linux', 'fold', 'all'): '6.97'}
ref_tup2sum3={(platform, 'fam', 'kappa'): 1.939, (platform, 'fam', 'all'): 1.94, (platform, 'sf', 'kappa'): 1.773, (platform, 'sf', 'all'): 1.806, (platform, 'fold', 'kappa'): 1.377, (platform, 'fold', 'all'): 1.436, ('linux', 'fam', 'kappa'): 1.939, ('linux', 'fam', 'all'): 1.94, ('linux', 'sf', 'kappa'): 1.773, ('linux', 'sf', 'all'): 1.806, ('linux', 'fold', 'kappa'): 1.376, ('linux', 'fold', 'all'): 1.436}
ref_tup2top3={(platform, 'sf', 'kappa'): 1.742, (platform, 'sf', 'all'): 1.747, (platform, 'fold', 'kappa'): 0.589, (platform, 'fold', 'all'): 0.586, ('linux', 'sf', 'kappa'): 1.742, ('linux', 'sf', 'all'): 1.747, ('linux', 'fold', 'kappa'): 0.589, ('linux', 'fold', 'all'): 0.586}

maxdt_pct = 10
maxmem_pct = 10
maxsum3_pct = 1
nerr = 0

def get_time_mem_from_log_file(fn):
	t = None
	m = None
	for line in open(fn):
		if line.startswith("Elapsed time "):
			t = line.split("Elapsed time ")[1].strip()
			smm, sss = t.split(':')
			t = int(smm)*60 + int(sss)
		if line.startswith("Max memory "):
			m = line.split("Max memory ")[1].strip().replace("Gb", "")
			m = float(m)
	return t, m

def get_sum3_from_log_file(fn):
	for line in open(fn):
		if line.find("Sum3=") > 0:
			match = re.search(r"Sum3=\s*(-?\d+\.\d+)", line)
			if match:
				return float(match.group(1))
			return None

def get_top3_from_log_file(fn):
	for line in open(fn):
		if line.find("Top3=") > 0:
			match = re.search(r"Top3=\s*(-?\d+\.\d+)", line)
			if match:
				return float(match.group(1))
			return None

tup2time = {}
tup2mem = {}
tup2sum3 = {}
tup2top3 = {}

# platforms = ("win", "linux")
platforms = (platform, )
truths = ("fam", "sf", "fold")
modes = ("kappa", "all")

if sfonly:
	truths = ("sf", )

for platform in platforms:
	for truth in truths:
		for mode in modes:
			# fold.all.win.search.log
			prefix = outdir + truth + "." + mode + "." + platform
			search_log_fn = prefix + ".search.log"
			sum3_log_fn = prefix + ".sum3.log"
			top3_log_fn = prefix + ".top3.log"

			t, m = get_time_mem_from_log_file(search_log_fn)
			sum3 = get_sum3_from_log_file(sum3_log_fn)
			top3 = None
			if truth != "fam":
				top3 = get_top3_from_log_file(top3_log_fn)

			tup = (platform, truth, mode)
			tup2time[tup] = t
			tup2mem[tup] = m
			tup2sum3[tup] = float(sum3)
			if top3:
				tup2top3[tup] = float(top3)

if 0:
	print(f"{tup2time=}")
	print(f"{tup2mem=}")
	print(f"{tup2sum3=}")
	print(f"{tup2top3=}")
	sys.exit(0)

for truth in truths:
	print()
	print("========= " + truth)
	if truth == "fam":
	  #     118( +0)  10.9(0.0)  1.94(+0.00)
	  print("  time      mem        sum3")
	else:
	  print("  time      mem        sum3         top3")
	for mode in modes:
		for platform in platforms:
			tup = (platform, truth, mode)
			t = tup2time[tup]
			ref_t = ref_tup2time[tup]
			dt = ref_t - t
			dtpct = dt*100/ref_t

			m = tup2mem[tup]
			ref_m = float(ref_tup2mem[tup])
			dm = ref_m - m
			mempct = dm*100/ref_m

			sum3 = tup2sum3[tup]
			ref_sum3 = ref_tup2sum3[tup]
			dsum3 = sum3 - ref_sum3
			sum3pct = -dsum3/ref_sum3

			s = "%5d" % t
			if dtpct > maxdt_pct:
				nerr += 1
				s += "(%+3d***)" % dt
			else:
				s += "(%+3d)" % dt

			s += "  %4.1f" % m
			if mempct > maxmem_pct:
				nerr += 1
				s += "(%3.1f***)" % dm
			else:
				s += "(%3.1f)" % dm

			s += "  %4.2f" % sum3
			if sum3pct > maxsum3_pct:
				nerr += 1
				s += "(%+3.2f***)" % dsum3
			else:
				s += "(%+3.2f)" % dsum3

			if truth != "fam":
				top3 = tup2top3.get(tup)
				ref_top3 = ref_tup2top3.get(tup)
				dtop3 = top3 - ref_top3
				s += "  %4.2f" % top3
				s += "(%+3.2f)" % dtop3

			print(s)
if nerr == 0:
	print("%s *** PASSED ***" % sys.argv[0])
	sys.exit(0)
else:
	print("%s *** FAILED ***" % sys.argv[0])
	sys.exit(1)
