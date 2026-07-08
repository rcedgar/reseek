def ReadEFAOnMSA(FileName, OnMSA):
	Labels = []
	Seqs = []
	Seq = ""
	Label = ""
	EFALabel = None
	for Line in open(FileName):
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == "<":
			if len(Seq) > 0:
				Seqs.append(Seq)
				Labels.append(Label)
			if len(Labels) > 0:
				OnMSA(EFALabel, Labels, Seqs)
			Labels = []
			Seqs = []
			EFALabel = Line[1:]
		else:
			if Line[0] == '>':
				if len(Seq) > 0:
					Labels.append(Label)
					Seqs.append(Seq)
				Label = Line[1:]
				Seq = ""
			else:
				Seq += Line.replace(" ", "")
	if len(Seq) > 0:
		Labels.append(Label)
		Seqs.append(Seq)
	if len(Labels) > 0:
		OnMSA(EFALabel, Labels, Seqs)

def ReadSeqsOnSeq(FileName, OnSeq):
	Label = None
	Seq = ""
	for Line in open(FileName):
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			if len(Seq) > 0:
				OnSeq(Label, Seq)
			Label = Line[1:]
			Seq = ""
		else:
			Seq += Line.replace(" ", "")
	if len(Seq) > 0:
		OnSeq(Label, Seq)

def ReadSeqsDict(FileName, trunclabels=False):
	SeqDict = {}
	Label = None
	Seq = ""
	for Line in open(FileName):
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			if len(Seq) > 0:
				SeqDict[Label] = Seq
			Label = Line[1:]
			if trunclabels:
				Label = Label.split()[0]
			Seq = ""
		else:
			Seq += Line.replace(" ", "")
	if len(Seq) > 0:
		SeqDict[Label] = Seq
	return SeqDict

def WriteSeq(File, Seq, Label, BLOCKLENGTH = 80):
	if len(Seq) == 0:
		return
	if Label != "":
		File.write(">" + Label + "\n")
	if BLOCKLENGTH <= 0:
		File.write(Seq + "\n")
	else:
		SeqLength = len(Seq)
		BlockCount = (SeqLength + (BLOCKLENGTH-1))//BLOCKLENGTH
		for BlockIndex in range(BlockCount):
			Block = Seq[BlockIndex*BLOCKLENGTH:]
			Block = Block[:BLOCKLENGTH]
			File.write(Block + "\n")

def GetAccFromLabel(Label):
	Fields = Label.split(';')
	if Fields[0] == "":
		return ""
	Fields2 = Fields[0].split()
	return Fields2[0]
