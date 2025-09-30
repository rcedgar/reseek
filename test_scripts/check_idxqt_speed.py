#!/usr/bin/python3

import socket
import sys
import os
import re

errors = 0

time_str = "idxq elapsedsecs=6.67 memkb=639796 pctcpu=88% systemsecs=0.25 usersecs=5.68"
time_regex = r'^([^ ]*) elapsedsecs=([0-9.]*) memkb=([0-9.]*) pctcpu=([0-9.]*)% systemsecs=([0-9.]*) usersecs=([0-9.]*)'

def do_time(fn):
	with open(fn) as f:
		line = f.readline()[:-1]
		M = re.match(time_regex, line)
		assert not M is None
		idx, secs, memkb, pctcpu, systemsecs, usersecs = M.groups()
		return float(secs), float(memkb)/1000

# 6.94 643264
# 86.23 639900

idxq_secs, idxq_Mb = do_time("../test_output/1hhs_pdb90_idxq.time")
idxt_secs, idxt_Mb = do_time("../test_output/1hhs_pdb90_idxt.time")
hostname = socket.gethostname()

if hostname == "rip":
	if idxq_secs > 8:
		errors += 1
		print("%s: ERROR 1hhs PDB90 idxq secs=%d (max 8)" \
			% (sys.argv[0], idxq_secs))
	if idxt_secs > 100:
		errors += 1
		print("%s: ERROR 1hhs PDB90 idxt secs=%d (max 100)" \
			% (sys.argv[0], idxt_secs))
	if idxq_Mb > 700:
		errors += 1
		print("%s: ERROR 1hhs PDB90 idxq Mb=%d (max 700)" \
			% (sys.argv[0], idxq_Mb))
	if idxt_Mb > 700:
		errors += 1
		print("%s: ERROR 1hhs PDB90 idxt Mb=%d (max 700)" \
			% (sys.argv[0], idxt_Mb))

if errors == 0:
	print("ok check_idxqt_speed")
exit(1 if errors > 0 else 0)
