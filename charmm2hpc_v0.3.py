# CHARMM2HPC.PY 
#this thing takes in a specified input file containing all necessary information
'''
The script has four modes:
taurus: 	BULLCLUSTER DD
hildi:		CLUSTER LE
umbrella:	umbrella mode for non-SLURM based run
default:	just correction of thermostats
'''
import os

version = "charmm2hpc.py v0.3 11/2019"
logfile = open("charmm2hpc.log","w+")
logfile.write("Initiating %s\n" % (version))
logfile.write("Correspondence: f.seufert@studserv.uni-leipzig.de\n\n")
activedir = os.getcwd()
logfile.write("Active Folder:"+activedir+"\n")
#LOAD DATA
'''
#INPUT FILE for charmm2hpc.py -> enter all parameters here >> charmm2hpc.inp
#Version 0.1 11/2019

simname		= 		bench_mult
duration	= 		5									#put your simduration in nanoseconds(1ns = 500,000 steps)
preset		= 		hildi								#put taurus for TAURUS(Dresden Bull cluster), put hildi for Leipzig cluster, umbrella for umbrella, anything else is unbiased normal sim
xtcsteps	= 		2500								#put interval of steps to be written to XTC-file
replicas	= 		3									#put number of replicas you want to generate, generates separate *.sh files for each run

#Input for PULL-Runs
biased		= 		N									#put yes or no, no skips the umbrella inclusion
pullgrp1	= 		Lauch1 44 162 388 589				#leave empty if no umbrella
pullgrp2	= 		Gemuese2 6001 6002 6027 6028		#leave empty if no umbrella
pullrate	= 		0.001
pullforce	= 		754
compressed	= 		Y									#put Y or N, umbrella output will then only be written to *.xtc file

#Input for BULLCLUSTER
username	= 		flse848b									
simtime		= 		96:00:00							#put hh;mm;ss, keep as low as necessary (taurus max 7 days)
mail		=		f.seufert@studserv.uni-leipzig.de
mailtype	= 		all									#all, error etc.
workdir		= 		/scratch/p_membranproteine/flo/gain			
node_number =		20									#20 should be sufficient for most systems <300k atoms
node_taskss =		16
nodetype	=		sandy								#put sandy or haswell for taurus
project		=		p_membranproteine

#Input for HILDI
username	=		fseufert									
interval/len=		2:00:00								#interval for staggered simulation
gpu			=		1									#leave at 1 ideally for GROMACS
tasks-p-node=		9                                   
cpus-p-task	=		2									#no. of cores per task
nice		=		0									#for running lower priority
mail-type	=		FAIL
'''

inputs=open("charmm2hpc.inp","r")
linect=0
input_data=[]
for line in inputs:
	if line.startswith("#")==True: continue
	if len(line) < 2: continue
	try:no_comment=line.split("#")
	except: no_comment=line
	content=no_comment[0].split("=")
	content[1] = content[1].strip(" \n\r\t")
#content[1] contains now the relevant information for each line
	input_data.append(content[1])
	content = []
logfile.write("\n Inputfile read.\n")
for item in input_data:
        logfile.write("%s\n" % item)
#transfer read data to variables
simname = input_data[0].rstrip(" \r\n\t")
try: duration = int(input_data[1])
except: logfile.write("ERROR:DURATION IS NOT AN INTEGER\n")
preset = input_data[2].rstrip(" \r\n")
if preset.startswith("taurus") == False and preset.startswith("hildi") == False: 
	logfile.write("NON-HPC PRESET: resorting to default mode\n")
try: xtcinterval = int(input_data[3])
except: logfile.write("ERROR:XTCSTEPS IS NOT AN INTEGER\n")
try: reps = int(input_data[4])
except: logfile.write("ERROR:REPLICAS IS NOT AN INTEGER\n")
if input_data[5].startswith("Y"): bias = True
elif input_data[5].startswith("N"): bias = False
else:
	bias = False
	logfile.write("Warning: Unrecognized character in Line 6. Assuming unbiased run\n")
if bias == True:	#read the pull block only when a "Y" is put into the biased row
	group1 = input_data[6].split(" ",1)
	group1_name = group1[0]
	group1_atoms = group1[1].rstrip(" \r\n")
	group2 = input_data[7].split(" ",1)
	group2_name = group2[0]
	group2_atoms = group2[1].rstrip(" \r\n")
	logfile.write("Biased run: Loaded umbrella data.\n"+group1_name+":"+group1_atoms+"\n"+group2_name+":"+group2_atoms)
	try:p_rate = float(input_data[8])
	except:logfile.write("ERROR:PULLRATE IS NOT FLOAT")
	try:p_force = int(input_data[9])
	except:logfile.write("ERROR:PULLRATE IS NOT INTEGER")
	if input_data[10].startswith("Y"): compressed_pull = True
	elif input_data[10].startswith("N"): compressed_pull = False
	else:
		compressed_pull = False
		logfile.write("Warning: Unrecognized character in Line 11. Assuming uncompressed biased run\n")

if preset.startswith("taurus") == True :
	bashdata = []
	for item in input_data[11:20]: bashdata.append(item)
	logfile.write("\nTAURUS run detected. Loading HPC section:\n")
	logfile.write(str(bashdata)+"\n\n")
if preset.startswith("hildi")==True:
	bashdata = []
	for item in input_data[20:28]: bashdata.append(item)
	logfile.write("\nHILDI run detected. Loading HPC section:\n")
	logfile.write(str(bashdata)+"\n\n")

#FUNCTION SECTION

def convert_linebreak(path):
	data = open(path, "rb").read()
	newdata = data.replace("\r\n", "\n")
	if newdata != data:
		logfile.write("Converted LineBreaks in BASH-File to UNIX\n")
		f = open("gromacs/%s.sh" % (simname), "wb")
		f.write(newdata)
		f.close()

def separator():
	logfile.write("\n-----------------------------------------------------------------------------\n\n")

def rewrite_step4_1(biased):
	separator()
	eq_mdp=open("gromacs\step4.1_equilibration.mdp", "r+")
	eq_new=open("gromacs\step4.1_%s.mdp" % (simname),"w+")
	logfile.write("Function rewrite_step4_1 called\n")
	replaced = 0
	line_nr = 0
	for line in eq_mdp:
		line_nr+=1
		if line_nr == 1 and biased == True: 
			eq_new.write(line.strip("\n"))
			eq_new.write(" -D%s -D%s\n" % (group1_name,group2_name))	#define groups for pulling
			logfile.write("defined %s and %s as pulling groups\n" % (group1_name,group2_name))
		elif line.startswith("tcoupl"): 
			eq_new.write(line[:26]+"Berendsen\n") #THIS NEEDS TO REPLACE THE WORD IN THE LINE!
			replaced = 1
		elif line_nr == 27 and biased == True:
			eq_new.write("pcoupl                  = Berendsen\npcoupltype              = semiisotropic\ntau_p                   = 5.0\ncompressibility         = 4.5e-5 4.5e-5\nref_p                   = 1.0    1.0\n;\n")
			logfile.write("Replaced pcoupl Block! Biased Sim!\n")
			eq_new.write(line)
		else: eq_new.write(line)
	if replaced == 0: logfile.write("\nWARNING: Could not find >Nose-Hoover< in step4.1_equilibration.mdp! Proceeding")
	else: logfile.write("\nReplaced tcoupl with >Berendsen< in step4.1_equilibration.mdp!\n")
	if biased == True:
		eq_new.write('''
gen-vel					= yes		; not necessarily needed in all multiplicates
gen-temp				= 310.15
gen-seed				= -1



; Periodic boundary conditions are on in all directions
pbc     = xyz
; Long-range dispersion correction
DispCorr    = EnerPres
; Pull code
pull                    = yes
pull_ncoords            = 1         ; only one reaction coordinate
pull_ngroups            = 2         ; two groups defining one reaction coordinate
pull_group1_name        = %s
pull_group2_name        = %s
pull_coord1_type        = umbrella  ; harmonic potential
pull_coord1_geometry    = distance  ; simple distance increase
pull_coord1_dim         = Y Y Y
pull_coord1_groups      = 1 2
pull_coord1_start       = yes       ; define initial COM distance > 0
pull_coord1_rate        = 0.000      ; in equilibration, keep distance constant
pull_coord1_k           = 1000      ; kJ mol^-1 nm^-2
pull-nstxout            = 10
pull-nstfout            = 1
;use -pf and -px names after the filename for data analysis
''' % (group1_name,group2_name))
		corrected_rate = str(p_rate * 100)
		logfile.write("Wrote PULL Block into file %s.mdp:\nG1: %s\nG2: %s\nRate: ZERO (equilibration) nm/nsec\nForce: 1000 kJ mol^-1 nm^-2\n" % (simname,group1_name,group2_name))

	eq_mdp.close()
	eq_new.close()
	
def rewrite_step5(biased=False):
	separator()
	logfile.write("Function rewrite_step5 called. Rewriting step5 based on BIAS = %s\n" % (str(biased)))
	prod_mdp=open("gromacs\step5_production.mdp", "r+")
	prod_new=open("gromacs\%s.mdp" % (simname), "w+")
	compressed = True
	if biased == True and compressed_pull == False: compressed = False#compression for Header: will only be written uncompressed if explicitly stated
	if biased == True:
		prod_new.write("define                  =  -D%s -D%s\n" % (group1_name,group2_name))	#define groups for pulling
		logfile.write("defined %s and %s as pulling groups\n" % (group1_name,group2_name))
	line_nr=0
	for line in prod_mdp:
		line_nr+=1
		if line_nr == 4: prod_new.write("nstlog                  = %d\n" % (xtcinterval))
		elif line_nr == 5 and compressed == True:
			prod_new.write('''nstxout                 = 0            ; Do not want .trr-files
nstvout                 = 0             ; No velocities in output
nstfout                 = 0            ; No forces in output\n''')
		elif line_nr > 5 and line_nr < 8 and compressed == True: logfile.write("Skipping line %d\n" % (line_nr))
		elif line_nr == 9:
			prod_new.write("nstenergy               = %d\n" % (xtcinterval))
			logfile.write("set up compressed interval: %d into XTC file\n" % (xtcinterval))
			prod_new.write("nstxout-compressed      = %d\n ; But do want .xtc-files" % (xtcinterval))       
		elif line.startswith("tcoupl"): 
			prod_new.write(line[:26]+"v-rescale\n")
		elif line.startswith("nsteps"):
			steps=500000*duration
			prod_new.write(line[:26]+"%d\n" % (steps))
			logfile.write("Simulation Duration will be %d ns. Total %d steps\n" % (duration,steps))
		elif line.startswith("continuation") and biased == True: continue #usually line 35
		elif line_nr == 26 and biased == True:
			prod_new.write("pcoupl                  = Berendsen\npcoupltype              = semiisotropic\ntau_p                   = 5.0\ncompressibility         = 4.5e-5 4.5e-5\nref_p                   = 1.0    1.0\n")
			logfile.write("Replaced pcoupl Block! Biased Sim!\n")
		elif line_nr > 26 and line_nr < 32 and biased == True: logfile.write("Skipping Line %d\n" % (line_nr))
		else: prod_new.write(line)
#From Here, insert the biased md parameters
	if biased == True:
		prod_new.write('''
gen-vel					= yes		; not necessarily needed in all multiplicates
gen-temp				= 310.15
gen-seed				= -1
;

; Periodic boundary conditions are on in all directions
pbc     = xyz
; Long-range dispersion correction
DispCorr    = EnerPres
; Pull code
pull                    = yes
pull_ncoords            = 1         ; only one reaction coordinate
pull_ngroups            = 2         ; two groups defining one reaction coordinate 
pull_group1_name        = %s
pull_group2_name        = %s
pull_coord1_type        = umbrella  ; harmonic potential
pull_coord1_geometry    = distance  ; simple distance increase 
pull_coord1_dim         = Y Y Y
pull_coord1_groups      = 1 2
pull_coord1_start       = yes       ; define initial COM distance > 0
pull_coord1_rate        = %s      ; 0.1 = 0.01 nm per ps = 10 nm per ns
pull_coord1_k           = %d      ; kJ mol^-1 nm^-2
pull-nstxout            = 10
pull-nstfout            = 1
;use -pf and -px names after the filename for data analysis
''' % (group1_name,group2_name,str(p_rate),p_force))
		corrected_rate = str(p_rate * 100)
		logfile.write("Wrote PULL Block into file %s.mdp:\nG1: %s\nG2: %s\nRate: %s nm per nsec\nForce: %d kJ mol^-1 nm^-2\n" % (simname,group1_name,group2_name,corrected_rate,p_force))
	prod_mdp.close()

def determine_indexchain(atomstring):
	pdbfile = open("gromacs\step3_charmm2gmx.pdb")
	atomls = []
	atomls = atomstring.split()
	atomlist = []
	for item in atomls:
		atomlist.append(int(item))
	logfile.write("determine_indexchain called: use atoms %s\n" % (atomstring))
	line_nr = 0
	flag = False
	for line in pdbfile:
		line_nr += 1
		if line_nr > 3:
			atomnr = int(line[6:11].lstrip())
			for atom in atomlist:
				if atomnr == atom:
					chain = line[72:76].rstrip()
					logfile.write("Found %s in charmm2gmx.pdb. Chain: %s\n" % (atom,chain))
					flag = True
		if flag == True: break	
	pdbfile.close()
	return chain
	
def rewrite_topol():	#For Pull runs, a new topology file will be generated that loads the corresponding *.itp files
	separator()
	logfile.write("Function rewrite_topol called: editing links to include *.itp files\n")
	top=open("gromacs/topol.top","r+")
	top_new=open("gromacs/topol_pull.top","w+")
	line_nr = 0
	for line in top:
		line_nr += 1
		if line.startswith("#") == True:
			top_new.write(line)
			if line.find(group1_chain) != -1:
				top_new.write("#ifdef %s\n#include \"toppar/%s.itp\" ; Pull-Group itp file in toppar/\n#endif\n" % (group1_name,group1_name))
			if line.find(group2_chain) != -1:
				top_new.write("#ifdef %s\n#include \"toppar/%s.itp\" ; Pull-Group itp file in toppar/\n#endif\n" % (group2_name,group2_name))
		elif line_nr == 10: top_new.write(line[:-1]+" Modified by %s\n" % (version))
		else: top_new.write(line)
	top.close()
	top_new.close()
	logfile.write("Successfully edited topol-file. New groups: %s %s , %s %s" % (group1_name, group1_chain, group2_name, group2_chain))
	newtop = os.path.join(activedir, "gromacs/topol_pull.top")
	deftop = os.path.join(activedir, "gromacs/topol.top")
	oldtop = os.path.join(activedir, "gromacs/topol_old.top")
	try: os.remove(oldtop)
	except: logfile.write("No preexisting topol_old.top found. Moving old version to topol_old.top\n")
	os.rename(deftop,oldtop)
	os.rename(newtop,deftop)
	
def rewrite_ndx():		#For Pull runs, generate a new NDX file with the two new groups
	separator()
	logfile.write("Function rewrite_ndx called: Adding new indexgroups into index_pull.ndx\n")
	ind = open("gromacs/index.ndx","r+")
	ind_new = open("gromacs/index_pull.ndx","w+")
	for line in ind: 
		if len(line)>1:ind_new.write(line)
	ind_new.write(" [%s]\n%s\n" % (group1_name,group1_atoms))
	ind_new.write(" [%s]\n%s\n" % (group2_name,group2_atoms))
	logfile.write("Editing successful\n")
	ind.close()
	ind_new.close()
	newndx = os.path.join(activedir, "gromacs/index_pull.ndx")
	defndx = os.path.join(activedir, "gromacs/index.ndx")
	oldndx = os.path.join(activedir, "gromacs/index_old.ndx")
	try: os.remove(oldndx)
	except: logfile.write("No preexisting index_old.ndx found. Moving old version to index_old.ndx\n")
	os.rename(defndx,oldndx)
	os.rename(newndx,defndx)

def generate_itp(name,atoms,chain):	#For Pull runs, generate a new ITP file for the particular group! Check the chain and adjust!
	separator()
	logfile.write("Function generate_itp called with arguments: %s, %s, " % (name, chain))
	logfile.write(atoms)
	logfile.write("\n")
	itp_new = open("gromacs/toppar/%s.itp" % (name),"w+")
	pdbfile = open("gromacs\step3_charmm2gmx.pdb")
	#write header first
	itp_new.write("; position restraints for PROT of Title\n \n [ position_restraints ]\n;  i funct       fcx        fcy        fcz\n")
	atlst = atoms.split()
	atomlist = []
	for item in atlst: atomlist.append(int(item))
	#correct the numbers in atomlist by subtracting the begin of the correspondent chain --> each chain starts with atom number 1
	line_nr = 0
	chains=dict()
	current_chain = ""
	for line in pdbfile:
		line_nr+=1
		if line_nr<=3: continue
		elif line_nr>100000: break
		elif line[72:76] == current_chain: continue
		else: 
			current_chain = line[72:76].rstrip()
			chains[current_chain]=line_nr-4
			logfile.write("Found chain: %s in line %d\n" % (current_chain,line_nr))
	logfile.write("The following chains were found:")
	logfile.write(str(chains))
	#now, correct the index number of atoms with by subtracting with the extracted integer!
	corrected_atoms = []
	for atom in atomlist:
		corrected_atoms.append(atom - chains[chain])
	logfile.write("\nThe corrected atom numbers are:")
	logfile.write(str(corrected_atoms))
	#write a new line for each atom in the corrected list
	for number in corrected_atoms:
		if len(str(number)) < 4:
			for i in range(4-len(str(number))): itp_new.write(r" ") #adjust the indentation for length
		itp_new.write("%d    1       1000       1000       1000\n" % number)
	itp_new.write("\n")
	logfile.write("\n%s.itp successfully generated in /gromacs/toppar!\n" % (name))
	itp_new.close()
	pdbfile.close()

def generate_bash(hpc_type,iterations):
	separator()
	logfile.write("Called Function: generate_bash with %s preset and %d iterations!\n" % (hpc_type, int(iterations)))
	bash_new=open("gromacs/%s.sh" % (simname), "w+")	
	if hpc_type.startswith("taurus")==True:	#add header section for TAURUS
		bash_new.write('''#!/bin/bash
#SBATCH -A %s #project
#SBATCH -J %s #user
#SBATCH -p %s # you can select between sandy and haswell "https://doc.zih.tu-dresden.de/hpc-wiki/bin/view/Compendium/SystemTaurus#Partitions"
#SBATCH --nodes=%s #
#SBATCH --ntasks-per-node=%s #tasks/ranks per node # It should be changed to 24 if you run on haswell
#SBATCH --cpus-per-task=1 #physical cores per all tasks/ranks
#SBATCH --mem=5000 #RAM in MB per node
#SBATCH --time=%s # maximum 7 days
#SBATCH --mail-type=%s
#SBATCH --mail-user=%s
#SBATCH --workdir=\"%s\"
#SBATCH --output=simulation-m%d.out
#SBATCH --error=simulation-m%d.err

module load modenv/scs5
module load impi/2018.4.274-iccifort-2019.1.144-GCC-8.2.0-2.31.1\n\n\n''' % (bashdata[8],bashdata[0],bashdata[7],bashdata[5],bashdata[6],bashdata[1],bashdata[3],bashdata[2],bashdata[4],iterations,iterations))
		tasks = int(bashdata[5])*int(bashdata[6])
		logfile.write("Total Number of tasks to launch on taurus: %d" % (tasks))
		bash_new.write('''IT=1
name=%s_$IT
gmx grompp -f step4.0_minimization.mdp -c step3_charmm2gmx.pdb -r step3_charmm2gmx.pdb -n index.ndx -o step4.0_$IT.tpr -p topol.top;
srun -n %d gmx_mpi mdrun -deffnm step4.0_$IT -dlb yes;
gmx grompp -f step4.1_%s.mdp -c step4.0_$IT.gro -r step4.0_$IT.gro -n index.ndx -o step4.1_$name.tpr -p topol.top;
srun -n %d gmx_mpi mdrun -deffnm step4.1_$name -dlb yes;
gmx grompp -f %s.mdp -c step4.1_$name.gro -r step4.1_$name.gro -o $name.tpr -n index.ndx -p topol.top;
srun -n %d gmx_mpi mdrun -deffnm $name -dlb yes;''' % (simname,tasks,simname,tasks,simname,tasks))

	elif hpc_type.startswith("hildi")==True:  #add header section for HILDI
		bash_new.write('''#!/bin/bash
#SBATCH --gres=gpu:%s            #specifiy the number of GPUs, for GROMACS leave at 1
#SBATCH --ntasks-per-node=%s     #number of tasks
#SBATCH --cpus-per-task=%s	     #number of cores per task: total number of cores = number of tasks x number of cores per task: for GROMACS 18 cores should be used per GPU
#SBATCH --mem=4096               #total memory the job can use, 4096 should be fine for most GROMACS jobs
                                 #if you want to use Rosetta this option is more useful: #SBATCH --mem-per-cpu={amount of memory in MB}
#SBATCH --mem-bind=local         #only use memory that is directly connected to the used cores
#SBATCH --nice=%s                #if you want to run a lot of jobs without interfering too much with jobs of other users use a value larger than 0, this way the jobs of other users will have higher priority than yours, the maximum value is 2147483645
#SBATCH --time=%s 		         #time limit of the job, format=hours:minutes:seconds
#SBATCH --job-name=%s
#SBATCH --mail-type=%s           #send email if job fails
		
''' % (bashdata[2],bashdata[3],bashdata[4],bashdata[5],bashdata[1],simname,bashdata[6]))
		timesplit = bashdata[1].split(":",1)
		time_hours=int(timesplit[0])
		logfile.write("Extracted time in hours: %d\n" % (time_hours))
		if time_hours < 1: 
			logfile.write("ERROR: PLEASE INCREASE SIMULATION TO AT LEAST 1 HOUR!\n")
			exit()
		time_init=str(float(time_hours)-0.2)
		time_sim=str(float(time_hours)-0.01)
		bash_new.write('''module load cuda/10.1.243
module load GROMACS/2019.4_GPU
IT=1

name=%s_$IT
scriptN=%s.sh
#simulation time
totalT=%d
#Time step
dT=0.002
#Target time (ns)
finalT=%d
#Time step into ns
dT=`echo "$dT * 1000" | bc`
''' % (simname,simname,time_hours,duration))
		bash_new.write(r'''dT=${dT%.*}
#Target number of steps (fs)
fSteps0=$((finalT*1000000))
fSteps=$((fSteps0/dT))
#Get actual time from last xtc file
numxtc=`ls -1 ${name}.xtc ${name}.part*.xtc 2>/dev/null | wc -l`
if [[ "$numxtc" -eq 0 ]]; then
 xtcS=${numxtc}
else
 xtcFile=`ls -1tr ${name}.xtc ${name}.part*.xtc 2>/dev/null | tail -1`
 if [[ "$numxtc" -eq 1 ]] && [[ "$xtcFile" = "${name}.xtc" ]]; then
  xtcS=${numxtc}
 else
  xtcS=`ls -1tr ${name}.xtc ${name}.part*.xtc 2>/dev/null | tail -1 | rev | cut -c5-8 | rev | sed 's/^0*//'`
 fi
fi
if [[ "$xtcS" -gt 0 ]]; then
 steps0=`gmx check -f ${xtcFile} 2>&1 | awk '/Precision/{getline; print}' | sed -n -e 's/^.*time //p'`
 steps0=${steps0%.*} #Remove decimal places
 steps0=$((steps0*1000)) #Into number of steps
 Ns=$((steps0/1000000)) #Into ns
 steps=$((steps0/dT)) #Account for time step
else
 steps=0
 Ns=0
fi

debugF=debug${xtcS}.log
echo "Last file: $xtcFile" > $debugF
echo "Final: $fSteps steps ($((fSteps0/1000000)) ns)" >> $debugF
echo "Current: $steps steps ($Ns ns)" >> $debugF

#If the final number of steps is lower than the actual time, start/keep simulating.
if [[ "$steps" -lt "$fSteps" ]]; then
 nsteps=$((fSteps - steps))
 nsteps0=$((fSteps0 - steps0))
 echo "I am launching $nsteps steps ($((nsteps0/1000000)) ns)" >> $debugF
 #Run
 if [[ "$xtcS" -eq "0" ]]; then
''')
		bash_new.write(
'''  gmx grompp -f step4.0_minimization.mdp -c step3_charmm2gmx.pdb -r step3_charmm2gmx.pdb -n index.ndx -o step4.0_$IT.tpr -p topol.top;
  gmx mdrun -nt 18 -deffnm step4.0_$IT -dlb yes;
  gmx grompp -f step4.1_%s.mdp -c step4.0_$IT.gro -r step4.0_$IT.gro -n index.ndx -o step4.1_$name.tpr -p topol.top;
  gmx mdrun -nt 18 -deffnm step4.1_$name -dlb yes;
  gmx grompp -f %s.mdp -c step4.1_$name.gro -r step4.1_$name.gro -o $name.tpr -n index.ndx -p topol.top;
  ''' % (simname,simname))
		bash_new.write(r"gmx mdrun -nt 18 -nb gpu -pme gpu -deffnm $name -dlb yes -cpo ${name}_$((xtcS + 1)) -maxh "+"%s -nsteps $nsteps -v >& log;" % (time_init)+"\n else\n  ")
		bash_new.write(r"gmx mdrun -nt 18 -nb gpu -pme gpu -deffnm $name -dlb yes -cpi ${name}_${xtcS} -cpo ${name}_$((xtcS + 1)) -noappend -maxh "+"%s -nsteps $nsteps -v >& log\n" % (time_sim))
		bash_new.write(''' fi
 #Re-launch script
 sleep 20
 sbatch $scriptN
else
 echo "Final number of steps reached ($steps)" >> $debugF
 exit
fi
''')
	else: 
		logfile.write("WARNING: Unrecognized HPC type. Generating bash script for running on local gmx machine")
		bash_new.write('''IT=1
name=%s_$IT
gmx mdrun -nt 4 -deffnm step4.0_$IT -dlb yes;
gmx grompp -f step4.1_%s.mdp -c step4.0_$IT.gro -r step4.0_$IT.gro -n index.ndx -o step4.1_$name.tpr -p topol.top;
gmx mdrun -nt 4 -deffnm step4.1_$name -dlb yes;
gmx grompp -f %s.mdp -c step4.1_$name.gro -r step4.1_$name.gro -o $name.tpr -n index.ndx -p topol.top;
gmx mdrun -nt 4 -deffnm $name -dlb yes -pin on;
''' % (simname,simname))
	bash_new.close()
	#convert line breaks!
	convert_linebreak("gromacs/%s.sh" % (simname))
	#section for iteration: Generate additional bash files with different simnames! Add an iterative number to them
	if iterations > 1:
		logfile.write("Iterations: %d #Generating iterative BASH files\n" % (iterations))
		#Find instances of SIMNAME in header, output section and append number
		#Add also iterative numbers for each run of the minimization and equilibration
		for i in range(1,iterations+1):
			logfile.write("Iteration %d: Creating new file %s_%d.sh\n" % (i,simname,i))
			bash_it=open("gromacs/%s_%d.sh" % (simname,i),"w+")
			bash_read=open("gromacs/%s.sh" % (simname),"r+")
			line_nr=0
			for line in bash_read:
				line_nr+=1
				if line_nr <= 14: 
					bash_it.write(line)
				elif line.startswith("IT")==True: bash_it.write("IT=%d\n" % (i))
				elif line.startswith("scriptN=")==True: bash_it.write("scriptN=%s_%d.sh\n" %(simname,i))
				else: bash_it.write(line)
			bash_it.close()
			convert_linebreak("gromacs/%s_%d.sh" % (simname,i))
			
rewrite_step4_1(bias)
rewrite_step5(bias)

if bias == True: 
	group1_chain=determine_indexchain(group1_atoms)
	group2_chain=determine_indexchain(group2_atoms)
	rewrite_topol()
	rewrite_ndx()
	generate_itp(group1_name,group1_atoms,group1_chain)
	generate_itp(group2_name,group2_atoms,group2_chain)
generate_bash(preset,reps)
logfile.close()

print("EXECUTION COMPLETE")