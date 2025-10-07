import RNA  # pyright: ignore[reportMissingImports]
import subprocess
import os
import uuid

def rnafold(rnaseq):
	if isinstance(rnaseq, str):
		fc  = RNA.fold_compound(rnaseq)
		return [fc.mfe()]
	elif isinstance(rnaseq, list):
		secondary_structures, energies = [], []
		for seq in rnaseq:
			fc  = RNA.fold_compound(seq)
			ss, e = fc.mfe()
			secondary_structures.append(ss)
			energies.append(e)
		return secondary_structures, energies
	else:
		print("The input should be eather string or list of strings")
		
def rna_seq_search(headers, sequences, tmp):

	# run cmscan and read search results
	#tmp = "tmp" #+ str(uuid.uuid4())
	os.makedirs(tmp, exist_ok=True)



	with open(f"{tmp}/query.fa", "w") as f:
		for header, sequence in zip(headers, sequences):
			f.write(f">{header}\n{sequence}\n")
	
	infernal_command = [
		'cmscan',
		'--cpu', '24',
		'--tblout', f'{tmp}/results.tbl',
		'--rfam', #'--cut_ga', 
		'rfam/Rfam.cm',
		f'{tmp}/query.fa'
	]

	cmscan_result = subprocess.run(infernal_command, check=True, capture_output=True, text=True)

	with open(f"{tmp}/results.cmscan", 'w') as f:
		f.write(cmscan_result.stdout)


	#read results into a dict  
	scores = {}
	rfams = {}
	target_names = {}
	evalues = {}

	for line in open(f'{tmp}/results.tbl'): 
		if line.startswith("#"):
			continue
		splited_line = line.split()
		score = float(splited_line[14])
		evalue = float(splited_line[15])
		id = splited_line[2]
		rfam = splited_line [1]
		target_name = splited_line[0]

		if id in evalues:
			if evalue < evalues[id]:
				scores[id] = score  
				rfams[id] = rfam
				target_names[id] = target_name
				evalues[id] = evalue
		else:
			scores[id] = score  
			rfams[id] = rfam
			target_names[id] = target_name
			evalues[id] = evalue


	# create a sorted list and add 0/None if a query does not have a hit
	score_list = []
	rfam_list = []
	target_name_list = []
	eval_list = []

	for header in headers:
		if header in scores:
			score_list.append(scores[header])
			eval_list.append(evalues[header])
			rfam_list.append(rfams[header])
			target_name_list.append(target_names[header])
		else:
			score_list.append(0)
			eval_list.append(0)
			rfam_list.append("-")
			target_name_list.append("-")

	#os.removedirs('tmp')
	return score_list, eval_list, rfam_list, target_name_list
