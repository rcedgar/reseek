# scop40 library functions

import sys

# in Gb
def get_memory_usage():
	import resource
	kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	return kb*1000/1e9

# Scop identifier looks like a.1.2.3 
#     a    1           2      3
# class.fold.superfamily.family
def getfold(scopid):
	flds = scopid.split('.')
	assert len(flds) == 4
	fold = flds[0] + "." + flds[1]
	return fold

def getsf(scopid):
	flds = scopid.split('.')
	assert len(flds) == 4
	sf = flds[0] + "." + flds[1] + "." + flds[2]
	return sf

class Scop40:
	def score1_is_better(self, score1, score2):
		assert not score1 is None and not score2 is None
		if self.se == 's':
			return score1 > score2
		elif self.se == 'e':
			return score1 < score2
		else:
			assert False

	def set_possible_tfs(self):
		self.nrdoms = len(self.doms)
		self.nrdompairs = self.nrdoms*self.nrdoms - self.nrdoms

		assert self.nrdoms == 11211
		assert self.nrdompairs == 125675310

		if self.level == "fam1":
			self.NT = 108718
			self.NI = 927820
			self.NF = 124638772
		elif self.level == "sf2":
			self.NT = 454766
			self.NI = 581772
			self.NF = 124638772
		elif self.level == "sf3":
			self.NT = 454766
			self.NI = 0
			self.NF = 125220544
		elif self.level == "sf4":
			self.NT = 346048
			self.NI = 690490
			self.NF = 124638772
		elif self.level == "fold5":
			self.NT = 581772
			self.NI = 454766
			self.NF = 124638772
		elif self.level == "fold6":
			self.NT = 1036538
			self.NI = 0
			self.NF = 124638772
		else:
			assert False, self.level

		assert self.NT + self.NI + self.NF == 125675310

	def read_dom2scopid(self, dom2scopid_fn):
		self.doms = set()
		self.fams = set()
		self.sfs = set()
		self.folds = set()
		self.fam2doms = {}
		self.sf2doms = {}
		self.fold2doms = {}
		self.dom2fam = {}
		self.dom2sf = {}
		self.dom2fold = {}

		for line in open(dom2scopid_fn):
			flds = line[:-1].split('\t')
			assert len(flds) == 2
			dom = flds[0]
			scopid = flds[1]
			assert dom not in self.doms
			self.doms.add(dom)

			fam = scopid
			sf = getsf(scopid)
			fold = getfold(scopid)

			self.dom2fam[dom] = scopid
			self.dom2sf[dom] = sf
			self.dom2fold[dom] = fold

			if not fam in self.fams:
				self.fam2doms[fam] = []
				self.fams.add(fam)
			self.fam2doms[fam].append(dom)

			if not sf in self.sfs:
				self.sf2doms[sf] = []
				self.sfs.add(sf)
			self.sf2doms[sf].append(dom)

			if not fold in self.folds:
				self.fold2doms[fold] = []
				self.folds.add(fold)
			self.fold2doms[fold].append(dom)

		self.fam2size = {}
		for fam in self.fams:
			self.fam2size[fam] = len(self.fam2doms[fam])
		self.sf2size = {}
		for sf in self.sfs:
			self.sf2size[sf] = len(self.sf2doms[sf])
		self.fold2size = {}
		for fold in self.folds:
			self.fold2size[fold] = len(self.fold2doms[fold])

	def eval_unsorted(self, qs, ts, plot_scores):
		if not self.quiet:
			sys.stderr.write("sorting...\n")
		nrhits = len(qs)
		assert nrhits > 10
		assert len(ts) == nrhits
		assert len(plot_scores) == nrhits
		v = [ (plot_scores[i], i) for i in range(nrhits) ]
		do_reverse = (self.se == "s")
		v_sorted = sorted(v, reverse=do_reverse)
		qs_sorted = []
		ts_sorted = []
		scores_sorted = []
		for _, i in v_sorted:
			q = qs[i]
			t = ts[i]
			if q == t:
				continue
			qs_sorted.append(q)
			ts_sorted.append(t)
			scores_sorted.append(plot_scores[i])
		if not self.quiet:
			sys.stderr.write("...done\n")
		self.eval_sorted(qs_sorted, ts_sorted, scores_sorted)

	# 0-based field nrs
	def read_file(self, fn, qfldnr, tfldnr, scorefldnr):
		self.qs = []
		self.ts = []
		self.scores = []
		if not self.quiet:
			sys.stderr.write("Reading %s...\n" % fn)
		for line in open(fn):
			flds = line[:-1].split('\t')
			q = flds[qfldnr]
			t = flds[tfldnr]
			if q == t:
				continue
			self.qs.append(q)
			self.ts.append(t)
			self.scores.append(float(flds[scorefldnr]))
		if not self.quiet:
			sys.stderr.write("...done\n")

	def eval_file(self, fn, qfldnr, tfldnr, scorefldnr, is_sorted):
		self.read_file(fn, qfldnr, tfldnr, scorefldnr)
		if is_sorted:
			self.eval_sorted(self.qs, self.ts, self.scores)
		else:
			self.eval_unsorted(self.qs, self.ts, self.scores)

	# 1=TP, 0=FP, -1=ignore
	def is_tp(self, q, t):
		q = q.split('/')[0]
		t = t.split('/')[0]
		if self.level == "fam1":
			qfam = self.dom2fam.get(q)
			tfam = self.dom2fam.get(t)
			if qfam == tfam:
				return 1
			qfold = self.dom2fold.get(q)
			tfold = self.dom2fold.get(t)
			if qfold != tfold:
				return 0
			return -1

		elif self.level == "sf2":
			qsf = self.dom2sf.get(q)
			tsf = self.dom2sf.get(t)
			if qsf is None or tsf is None:
				return -1
			if qsf == tsf:
				return 1
			else:
				return 0

		elif self.level == "sf3":
			qsf = self.dom2sf.get(q)
			tsf = self.dom2sf.get(t)
			if qsf == tsf:
				return 1
			qfold = self.dom2fold.get(q)
			tfold = self.dom2fold.get(t)
			if qfold != tfold:
				return 0
			return -1

		elif self.level == "sf4":
			qfam = self.dom2fam.get(q)
			tfam = self.dom2fam.get(t)
			if qfam == tfam:
				return -1
			qsf = self.dom2sf.get(q)
			tsf = self.dom2sf.get(t)
			if qsf == tsf:
				return 1
			qfold = self.dom2fold.get(q)
			tfold = self.dom2fold.get(t)
			if qfold != tfold:
				return 0
			return -1

		elif self.level == "fold5":
			qsf = self.dom2sf.get(q)
			tsf = self.dom2sf.get(t)
			if qsf == tsf:
				return -1
			qfold = self.dom2fold.get(q)
			tfold = self.dom2fold.get(t)
			if qfold == tfold:
				return 1
			else:
				return 0

		elif self.level == "fold6":
			qfold = self.dom2fold.get(q)
			tfold = self.dom2fold.get(t)
			if qfold == tfold:
				return 1
			else:
				return 0

		assert False

	def eval_sorted(self, qs, ts, scores):
		self.qs = qs
		self.ts = ts
		self.scores = scores
		unknown_doms = set()
		nrhits = len(qs)
		assert nrhits > 10
		assert len(ts) == nrhits
		assert len(scores) == nrhits
		last_score = None

		self.ntp = 0         # accumulated nr of TP hits
		self.nfp = 0         # accumulated nr of FP hits

		self.ntpe1 = 0       # accumulated nr of TP hits with E<=1
		self.nfpe1 = 0       # accumulated nr of FP hits with E<=1

		tpstep = 0.01   # bin size for TPs (X axis tick marks)
		tprt = 0.01     # current TPR threshold, +=tpstep during scan

		self.plot_tprs = []
		self.plot_fprs = []
		self.plot_epqs = []
		self.plot_scores = []
		self.plot_precisions = []

		self.tpr_at_fpepq0_1 = None
		self.tpr_at_fpepq1 = None
		self.tpr_at_fpepq10 = None

		self.dom2score_firstfp = {}
		self.dom2score_firsttp = {}
		for dom in self.doms:
			self.dom2score_firstfp[dom] = None
			self.dom2score_firsttp[dom] = None

		self.tps = []
		ni = 0
		if not self.quiet:
			sys.stderr.write("scanning hits...\n")
		for i in range(nrhits):
			if not self.quiet:
				if i%100000 == 0:
					gb = get_memory_usage()
					pct = 100*i/nrhits
					sys.stderr.write("%.1f%%, %.1f Gb RAM used   \r" % (pct, gb))
			q = qs[i].split('/')[0]
			t = ts[i].split('/')[0]
			if q == t:
				assert False, "Self-hits not removed"
			score = scores[i]

			if i > 0 and last_score != score:
				if not self.score1_is_better(last_score, score):
					assert False, \
						f"Not sorted correctly {q=} {t=} {self.se=} {last_score=} {score=}"
			last_score = score

			tp = self.is_tp(q, t)
			if tp == 1:
				self.ntp += 1
				if self.se == 'e' and score <= 1:
					self.ntpe1 += 1
				if self.dom2score_firsttp.get(q) is None or \
					self.score1_is_better(score, self.dom2score_firsttp[q]):
					self.dom2score_firsttp[q] = score
			elif tp == 0:
				self.nfp += 1
				if self.se == 'e' and score <= 1:
					self.nfpe1 += 1
				if self.dom2score_firstfp.get(q) is None or \
					self.score1_is_better(score, self.dom2score_firstfp[q]):
					self.dom2score_firstfp[q] = score
			elif tp == -1:
				ni += 1
			else:
				assert False
			self.tps.append(tp)

			# tpr=true-positive rate
			tpr = float(self.ntp)/self.NT

			# fpepq = false-positive errors per query
			fpepq = float(self.nfp)/self.nrdoms

			# precision = tp/(tp + fp)
			precision = 0
			fpr = 0
			if self.ntp+self.nfp > 0:
				precision = self.ntp/(self.ntp + self.nfp)
				fpr = self.nfp/(self.ntp + self.nfp)

			if fpepq >= 0.1 and self.tpr_at_fpepq0_1 is None:
				self.tpr_at_fpepq0_1 = tpr
			if fpepq >= 1 and self.tpr_at_fpepq1 is None:
				self.tpr_at_fpepq1 = tpr
			if fpepq >= 10 and self.tpr_at_fpepq10 is None:
				self.tpr_at_fpepq10 = tpr
			if tpr >= tprt:
				self.plot_tprs.append(tprt)
				self.plot_fprs.append(fpr)
				self.plot_epqs.append(fpepq)
				self.plot_scores.append(last_score)
				self.plot_precisions.append(precision)

				tprt += tpstep

			last_score = score
		sys.stderr.write("ni=%d\n" % ni)
		if self.tpr_at_fpepq0_1 is None:
			self.tpr_at_fpepq0_1 = tpr
		if self.tpr_at_fpepq1 is None:
			self.tpr_at_fpepq1 = tpr
		if self.tpr_at_fpepq10 is None:
			self.tpr_at_fpepq10 = tpr

		nrhits = len(qs)
		assert len(self.tps) == nrhits

		if not self.quiet:
			sys.stderr.write("%.1f M hits, %.1f Gb RAM used   \n" % (nrhits/1e6, get_memory_usage()))

		if self.tpr_at_fpepq0_1 is None:
			self.tpr_at_fpepq0_1 = tpr

		if self.tpr_at_fpepq1 is None:
			self.tpr_at_fpepq1 = tpr

		if self.tpr_at_fpepq10 is None:
			self.tpr_at_fpepq1 = tpr
	
		tpr = float(self.ntp)/self.NT
		fpepq = float(self.nfp)/self.nrdoms

		self.plot_tprs.append(tprt)
		self.plot_fprs.append(fpr)
		self.plot_epqs.append(fpepq)
		self.plot_scores.append(last_score)
		self.plot_precisions.append(precision)

		self.nrtps_to_firstfp = 0
		for i in range(nrhits):
			score = scores[i]
			tp = self.tps[i]
			q = qs[i].split('/')[0]
			tp = self.tps[i]
			if tp and (self.dom2score_firstfp.get(q) is None or \
				self.score1_is_better(score, self.dom2score_firstfp[q])):
				self.nrtps_to_firstfp += 1
		self.sens_to_firstfp = float(self.nrtps_to_firstfp)/self.NT
		nr_unknown = len(unknown_doms)
		if nr_unknown > 0:
			sys.stderr.write("Warning %d unknown doms\n" % nr_unknown)
			s = ""
			k = 0
			for dom in unknown_doms:
				s += " " + dom
				k += 1
				if nr_unknown > k:
					s += "..."
					break
			sys.stderr.write(s + "\n")

	def __init__(self, se, level, dom2scopid_fn, quiet = False):
		assert se == "s" or se == "e" # score or E-value
		self.se = se
		self.level = level
		self.quiet = quiet
		if se == 's':
			self.low_score = -1
		elif se == 'e':
			self.low_score = 99999
		else:
			assert False

		self.unknown_doms = set()
		self.read_dom2scopid(dom2scopid_fn)
		self.set_possible_tfs()

	def plot2file(self, fn):
		f = open(fn, "w")
		self.roc2filehandle(f)
		f.close()

	def plot2filehandle(self, f):
		n = len(self.plot_tprs)
		if len(self.plot_fprs) != n:
			assert False, "%d,%d" % (len(self.plot_fprs), n)
		assert len(self.plot_epqs) == n
		assert len(self.plot_scores) == n
		f.write("tpr\tepq\tfpr\tprecision\tscore\n")
		for i in range(n):
			s = "%.4g" % self.plot_tprs[i]
			s += "\t%.4g" % self.plot_epqs[i]
			s += "\t%.4g" % self.plot_fprs[i]
			s += "\t%.4g" % self.plot_precisions[i]
			s += "\t%.4g" % self.plot_scores[i]
			f.write(s + '\n')

	def roc_area(self, lo_epq, hi_epq):
		n = len(self.plot_tprs)
		assert len(self.plot_epqs) == n
		assert len(self.plot_scores) == n
		total = 0
		for i in range(n):
			epq = self.plot_epqs[i]
			if epq >= lo_epq and epq <= hi_epq:
				tpr = self.plot_tprs[i]
				total += tpr
		return total

	def get_summary(self):
		area = self.roc_area(0.01, 10)
		summary = "SEPQ0.1=%.4f" % self.tpr_at_fpepq0_1
		summary += " SEPQ1=%.4f" % self.tpr_at_fpepq1
		summary += " SEPQ10=%.4f" % self.tpr_at_fpepq10
		summary += " S1FP=%.4f" % self.sens_to_firstfp
		summary += " N1FP=%u" % self.nrtps_to_firstfp
		summary += " area=%.3g" % area
		return summary

	def get_summary2(self):
		summary = "TP=%u" % self.ntp
		summary += " TP1=%u" % self.ntpe1
		summary += " FP=%u" % self.nfp
		summary += " FP1=%u" % self.nfpe1
		return summary

	def top_hit_report1(self, f, dom, scoretp, scorefp, ntp, nfp):
		nr_doms = len(self.doms)
		nr_singleton_sfs = 838
		nr_possible_tps = nr_doms - nr_singleton_sfs
		s = dom
		if scoretp is None:
			s += "\t."
		else:
			s += "\t%.3g" % scoretp

		if scorefp is None:
			s += "\t."
		else:
			s += "\t%.3g" % scorefp
		s += "\t%.4f\t%.4f" % (ntp/nr_possible_tps, nfp/nr_doms)
		f.write(s + "\n")

	def top_hit_report(self, fn):
		K = 100
		f = open(fn, "w")
		ntp = 0
		nfp = 0
		v = []
		for dom in self.doms:
			scoretp = self.dom2score_firsttp[dom]
			scorefp = self.dom2score_firstfp[dom]
			if scoretp is None and scorefp is None:
				continue
			if scoretp is None and scorefp is None:
				if self.score1_is_better(scoretp, scorefp):
					sort_score = scoretp
				else:
					sort_score = scorefp
			elif scoretp is None:
				assert not scorefp is None
				sort_score = scorefp
			else:
				assert not scoretp is None
				sort_score = scoretp
			pair = (sort_score, dom)
			v.append(pair)
		do_reverse = (self.se == "s")
		v_sorted = sorted(v, reverse=do_reverse)
		for score, dom in v_sorted:
			scoretp = self.dom2score_firsttp[dom]
			scorefp = self.dom2score_firstfp[dom]
			if scoretp is None and scorefp is None:
				continue
			elif scoretp is None and not scorefp is None:
				nfp += 1
			elif not scoretp is None and scorefp is None:
				ntp += 1
			elif not scoretp is None and not scorefp is None:
				if self.score1_is_better(scoretp, scorefp):
					ntp += 1
				else:
					nfp += 1
			else:
				assert False
			if (ntp + nfp)%K == 0:
				self.top_hit_report1(f, dom, scoretp, scorefp, ntp, nfp)
		self.top_hit_report1(f, dom, scoretp, scorefp, ntp, nfp)
		f.close()

	def report_standard_counts(self):
		for self.level in [ "fam1", "sf2", "sf3", "sf4", "fold5", "fold6" ]:
			NT = 0
			NI = 0
			NF = 0
			k = 0
			for dom1 in self.doms:
				k += 1
				if k%100 == 0:
					pct = k*100.0/self.nrdoms
					sys.stderr.write("%.1f%%\r" % pct)
				for dom2 in self.doms:
					if dom1 == dom2:
						continue
					n = self.is_tp(dom1, dom2)
					if n == 1:
						NT += 1
					elif n == 0:
						NF += 1
					elif n == -1:
						NI += 1
					else:
						assert False
			print(f"{self.level=}  {NT=}  {NI=}  {NF=}  {NT+NI+NF=}")
			# NT1=  108718  NI1=  927820  NF1=124638772  N1=125675310
			# NT2=  454766  NI2=       0  NF2=125220544  N2=125675310
			# NT3=  454766  NI3=  581772  NF3=124638772  N3=125675310
			# NT4=  346048  NI4=  690490  NF4=124638772  N4=125675310
			# NT5=  581772  NI5=  454766  NF5=124638772  N5=125675310
			# NT6= 1036538  NI6=       0  NF6=124638772  N6=125675310
