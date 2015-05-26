#!/bin/bash

# Script for automating stochastic simulations of zebrafish segmentation
# Copyright (C) 2012 Ahmet Ay, Jack Holland, Adriana Sperlea

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# $1 = a file with parameter sets separated by line
# $2 = the directory where the data is or will be (if -J then this directory must already contain the simulation concentrations)
# -S, --scratch = use scratch space
# -J, --just-analyze = just analyze; do not run the simulations (the simulations must already exist)
# -F, --find-parameters = find parameters (i.e. stop if a set does not meet the proper criteria)
# -D, --delete-failed = delete directories for parameters that fail (only applies when -F is used)
# -q, --quiet = hide terminal output
# $...$n = the list of arguments to be passed to the simulation
#          (if -r|--runs or -m|--minutes are given as arguments to the simulation then this script uses them for analysis, otherwise they default to 1 and 600, respectively)

# this function prints out how to use this script and then exits
function usage () {
	echo "Usage: seg-clock <input file> <output directory> [-option[option]... [value]]... [--option [value]]..."
	exit 0
}

# this function prints out how to use this script and then exits
function echo () {
	if [ $quiet -eq 0 ]; then
		command echo "$1" "$2"
	fi
}

# this function submits a simulation instance to the job queue
function run_sim () {
	mkdir -p $rootdir/par$par
	
	echo "
	#PBS -N stochastic-$par
	#PBS -l nodes=1:ppn=1
	#PBS -l mem=1gb
	#PBS -l file=1gb
	#PBS -l walltime=48:00:00
	#PBS -q biomath
	#PBS -m n
	#PBS -j oe
	#PBS -o $rootdir/par$par/terminal-output.txt
	
	function check_dir_exists () {
	# \$1 mutant type directory name
		
		if [ ! -d \$pardir/\$1 ]; then
			if [ $findpars -eq 1 ]; then
				echo \"\$pardir/\$1 does not exist. ${color_red}Parameter set deemed invalid. Parameter set simulations stopped.${color_reset}\"
			else
				echo \"No data for \$1\"
			fi
			save_results
			if [ $findpars -eq 1 ]; then
				exit 1
			else
				exit 0
			fi
		fi
	}
	
	function echo () {
	# \$1 command to echo or option to echo
	# \$2 command to echo if \$1 was an option
		
		if [ $quiet -eq 0 ]; then
			command echo \"\$1\" \"\$2\"
		fi
	}
	
	function ofeatures () {
	# \$1 mutant type directory name
	# \$2-\$3 acceptable range for the mutant period
	# \$4 mutant type full name
		
		echo -n \"${color_blue}Analyzing oscillation features ${color_reset}... \"
		period_avg=0.0
		amp_avg=0.0
		ptot_avg=0.0
		
		check_dir_exists \$1
		
		for (( i=0; i < $runs; i++ ))
		do		
			values=\`\$curdir/analysis/ofeatures \$pardir/\$1/run\${i}/run\${i}.txt \$pardir/\$1/run\${i}/run\${i}_smooth.txt 100\`
			if [[ \"\$values\" = *[:digit:]* ]]; then
				behavior=\"\$behavior,error\"
				save_results
				exit 1
			fi
			echo \$values
			period=\${values%%,*}
			values=\${values#*,}
			amplitude=\${values%%,*}
			values=\${values#*,}
			ptot=\${values%%,*}

			period_avg=\$(\$curdir/analysis/calc \"\$period + \$period_avg\")
			amp_avg=\$(\$curdir/analysis/calc \"\$amplitude + \$amp_avg\")
			ptot_avg=\$(\$curdir/analysis/calc \"\$ptot + \$ptot_avg\")
		done
		
		period_avg=\$(\$curdir/analysis/calc \"\$period_avg / $runs\")
		amp_avg=\$(\$curdir/analysis/calc \"\$amp_avg / $runs\")
		ptot_avg=\$(\$curdir/analysis/calc \"\$ptot_avg / $runs\")
		
		behavior=\"\$behavior,\$period_avg,\$amp_avg\"
		if [ \$1 == \"wt\" ]; then
			period_wt=\$period_avg
			echo -n $output_done
		elif [ $findpars -eq 1 ]; then
			ratio=\$(\$curdir/analysis/calc \"\$period_avg / \$period_wt\")
			if [ \$(\$curdir/analysis/calc \"\$ratio < \$2\") -eq 1 ] || [ \$(\$curdir/analysis/calc \"\$ratio > \$3\") -eq 1 ]; then
				echo -n \"The \$4 vs wild type period ratio is not between \$2 and \$3. ${color_red}Parameter set deemed invalid. Parameter set simulations stopped.${color_reset}\"
				behavior=\"\$behavior,failed\"
				rm -rf \$pardir
				save_results
				exit 1
			else
				echo -n \"${color_blue}OK${color_reset}: The \$4 vs wild type period ratio is between \$2 and \$3\"
			fi
		fi
		echo \" (period=\$period_avg, wild type period=\$period_wt, amplitude=\$amp_avg, peak to trough ratio=\$ptot_avg)\"
	}
	
	function run_mut () {
	# \$1 = mutant type directory name
	# \$2 = mutant type full name
		
		if [ $justanalyze -eq 0 ]; then
			echo -n \"${color_blue}Creating ${color_reset}\$pardir/\$1 ... \"
			mkdir -p \$pardir/\$1
			echo $output_done
			echo -n \"${color_blue}Creating ${color_reset}\$pardir/\$1/par$par.txt ... \"
			if [ \$1 == \"wt\" ]; then
				pars=\"$1,$2,$3,$4,$5\"
			elif [ \$1 == \"her1\" ]; then
				pars=\"0.0,$2,$3,$4,$5\"
			elif [ \$1 == \"her7\" ]; then
				pars=\"$1,0.0,$3,$4,$5\"
			elif [ \$1 == \"her13\" ]; then
				pars=\"$1,$2,0.0,$4,$5\"
			elif [ \$1 == \"her713\" ]; then
				pars=\"$1,0.0,0.0,$4,$5\"
			else
				pars=\"$1,$2,$3,0.0,$5\"
			fi
			echo \$pars > \$pardir/\$1/par$par.txt
			echo $output_done
			
			echo \"${color_blue}Running ${color_reset}\$2 simulation ...\"
			time \$curdir/stochastic $args -i \$pardir/\$1/par$par.txt -o \$pardir/\$1
			if [ $? -ne 0 ]; then
				behavior=\"\$behavior,error\"
				save_results
				exit 1
			fi

			for (( i=0; i<$runs; i++ ))
			do
				echo -n \"${color_blue}Creating ${color_reset}\$pardir/\$1/run\${i} ... \"
				mkdir -p \$pardir/\$1/run\${i}
				echo $output_done
				echo -n \"${color_blue}Moving ${color_reset}\$pardir/\$1/run\${i}.txt to \$pardir/\$1/run\${i}/run\${i}.txt ... \"
				mv \$pardir/\$1/run\${i}.txt \$pardir/\$1/run\${i}/
				echo $output_done
			done

			echo \"${color_blue}Done simulating ${color_reset}\$2\"
			echo -n \"${color_blue}Smoothing data ${color_reset}... \"
			 
			for (( i=0; i < $runs; i++ ))
			do
				\$curdir/analysis/smoothing \$pardir/\$1/run\${i}/run\${i}.txt \$pardir/\$1/run\${i}/run\${i}_smooth.txt 40
			done		
			echo $output_done
		fi
	}
	
	function save_results () {
		echo -n \"${color_blue}Writing to ${color_reset}\$curdir/$rootdir/behavior.csv ... \"
		(
			flock -x 200
			echo \"\$behavior\">>\$curdir/$rootdir/behavior.csv
		) 200>>/dev/null
		echo $output_done
		
		if [ $delete_scratch -eq 1 ]; then
			echo -n \"${color_blue}Moving scratch directory ${color_reset}(\$pardir) to output directory (\$curdir/$rootdir) ... \"
			cp -r \$pardir \$curdir/$rootdir
			rm -rf \$pardir
			echo $output_done
		fi
		
		time_end=\$(date +%s)
		echo -e \"\nTotal elapsed time: \$((\$time_end - \$time_start))s\"
	}

	function sync() {
	# \$1 = mutant type
		
		total=0.0
		for (( i=0; i<$runs; i++ ))
		do
			syncscore=\`python \$curdir/analysis/synchronized.py \$pardir/\$1/run\${i}/run\${i}.txt $minsync\`
			if [ \${#syncscore} -gt 13 ]; then
				behavior=\"\$behavior,error\"
				save_results
				exit 1
			fi
			total=\$(\$curdir/analysis/calc \"\$syncscore + \$total\")
		done
		total=\$(\$curdir/analysis/calc \"\$total / $runs\")
		if [ \$1 == \"wt\" ]; then
			syncscore_wt=\$total
		elif [ \$1 == \"delta\" ]; then
			syncscore_delta=\$total
		elif [ \$1 == \"her1\" ]; then
			syncscore_her1=\$total
		elif [ \$1 == \"her7\" ]; then
			syncscore_her7=\$total
		elif [ \$1 == \"her13\" ]; then
			syncscore_her13=\$total
		else
			syncscore_her713=\$total
		fi
	}
	
	time_start=\$(date +%s)
	trap save_results SIGTERM
	curdir=\$PBS_O_WORKDIR
	pardir=$scratch/$rootdir/par$par
	behavior=\"$par\"
	if [ $justanalyze -eq 0 ]; then
		echo -n \"${color_blue}Creating ${color_reset}\$pardir if necessary ... \"
		mkdir -p \$pardir
		echo $output_done
	else
		echo -n \"${color_blue}Checking if ${color_reset}\$pardir exists ... \"
		if [ -d \$pardir ]; then
			echo $output_done
		else
			echo \"${color_red}Not found!${color_reset}\"
			behavior=\"\$behavior,error\"
			save_results
			exit 1
		fi
	fi
	
	echo \"${color_blue}Starting ${color_reset}wild type ... \"
	run_mut wt \"wild type\"
	syncscore_wt=0.0
	syncscore_delta=0.0
	syncscore_her1=0.0
	syncscore_her7=0.0
	syncscore_her13=0.0
	syncscore_her713=0.0
	echo -n \"${color_blue}Evaluating synchronization between cells ${color_reset}... \"
	sync wt
	behavior=\"\$behavior,\$syncscore_wt\"
	if [ $findpars -eq 1 ] && [ \$(\$curdir/analysis/calc \"\$syncscore_wt > 0.8\") -eq 0 ]; then
		echo \"Cells are not synchronized in the wild type (syncscore=\$syncscore_wt). ${color_red}Parameter set deemed invalid. Parameter set simulations stopped.${color_reset}\"
		behavior=\"\$behavior,failed\"
		rm -rf \$pardir
		save_results
		exit 1
	else
		if [ $findpars -eq 1 ]; then
			echo \"${color_blue}OK${color_reset}: Cells are synchronized in the wild type (syncscore=\$syncscore_wt)\"
		fi
		ofeatures wt 0 0 \"wild type\"
	fi
	echo \"${color_blue}Done with ${color_reset}wild type\"
	
	echo \"${color_blue}Starting ${color_reset}delta mutant ... \"
	run_mut delta \"delta mutant\"
	echo -n \"${color_blue}Evaluating synchronization between cells ${color_reset}... \"
	sync delta		
	behavior=\"\$behavior,\$syncscore_delta\"
	if [ $findpars -eq 1 ] && [ \$(\$curdir/analysis/calc \"\$syncscore_delta < 0.7\") -eq 0 ]; then
		echo \"Cells do not fall out of sync in the delta mutant (syncscore=\$syncscore_delta). ${color_red}Parameter set deemed invalid. Parameter set simulations stopped.${color_reset}\"
		behavior=\"\$behavior,failed\"
		rm -rf \$pardir
		save_results
		exit 1
	else
		if [ $findpars -eq 1 ]; then
			echo \"${color_blue}OK${color_reset}: Cells fall out of sync in the delta mutant (syncscore=\$syncscore_delta)\"
		fi
		ofeatures delta 1.07 1.30 \"Delta mutant\"
	fi
	echo \"${color_blue}Done with ${color_reset}Delta mutant\"
	
	echo \"${color_blue}Starting ${color_reset}Her13 mutant ... \"
	run_mut her13 \"Her13 mutant\"
	ofeatures her13 1.03 1.09 \"Her13 mutant\"
	echo \"${color_blue}Done with ${color_reset}Her13 mutant\"
	
	echo \"${color_blue}Starting ${color_reset}Her1 mutant ... \"
	run_mut her1 \"Her1 mutant\"
	ofeatures her1 0.97 1.03 \"Her1 mutant\"
	echo \"${color_blue}Done with ${color_reset}Her1 mutant\"
	
	echo \"${color_blue}Starting ${color_reset}Her7 mutant ... \"
	run_mut her7 \"Her7 mutant\"
	ofeatures her7 0.97 1.03 \"Her7 mutant\"
	echo \"${color_blue}Done with ${color_reset}Her7 mutant\"
	
	echo \"${color_blue}Starting ${color_reset}Her713 mutant ... \"
	run_mut her713 \"Her713 mutant\"
	ofeatures her713 1.03 1.09 \"Her713 mutant\"
	echo \"${color_blue}Done with ${color_reset}Her713 mutant\"
	
	if [ $findpars -eq 1 ]; then
		echo \"Wild type and all mutant conditions satisfied. ${color_green}Parameter set deemed valid.${color_reset}\"
	fi
	echo -n \"${color_blue}Calculating synchronization between cells ${color_reset} for Her1, Her7, Her13, and Her713... \"
	sync her1
	sync her7
	sync her13
	sync her713
	behavior=\"\$behavior,\$syncscore_her1,\$syncscore_her7,\$syncscore_her13,\$syncscore_her713\"
	echo $output_done
	
	if [ $findpars -eq 1 ]; then
		echo -n \"${color_blue}Writing to ${color_reset}\$curdir/$rootdir/allpassed.csv ... \"
		(
			flock -x 200
			echo \"$par,$1,$2,$3,$4,$5\">>\$curdir/$rootdir/allpassed.csv
		) 200>>/dev/null
		echo $output_done
	fi
	
	behavior=\"\$behavior,passed\"
	save_results
	" > "$rootdir/par$par/pbs-job"
	qsub "$rootdir/par$par/pbs-job"
}

# option to quiet terminal output
quiet=0

# output colors and shortcuts
color_blue=$(tput setaf 4)
color_green=$(tput setaf 2)
color_red=$(tput setaf 1)
color_reset=$(tput sgr0)
output_done=\"${color_blue}Done${color_reset}\"

# make sure the user provided proper arguments
if [ $# -lt 2 ] || [ ! -e "$1" ]; then
	usage
fi

parsets=$1
rootdir=$2
mkdir -p $rootdir

# create a file to store the parameter sets that passed
echo "set,psh1,psh7,psh13,psd,pdh1,pdh7,pdh13,pdd,msh1,msh7,msh13,msd,mdh1,mdh7,mdh13,mdd,ddgh1h1,ddgh1h7,ddgh1h13,ddgh7h7,ddgh7h13,ddgh13h13,delaymh1,delaymh7,delaymh13,delaymd,delayph1,delayph7,delayph13,delaypd,dah1h1,ddh1h1,dah1h7,ddh1h7,dah1h13,ddh1h13,dah7h7,ddh7h7,dah7h13,ddh7h13,dah13h13,ddh13h13,critph1h1,critph7h13,critpd" > $rootdir/allpassed.csv
echo "set,syncscore wt,per wt,amp wt,syncscore delta,per delta,amp wt,per her13,amp her13,per her1,amp her1,per her7,amp her7,per her713,amp her713,syncscore her1,syncscore her7,syncscore her13,syncscore her713" > $rootdir/behavior.csv

# check and process user arguments into useful variables
scratch="\$PBS_O_WORKDIR"
delete_scratch=0
justanalyze=0
findpars=0
deletefailed=0
runs=1
minsync=600
args=""

lastok=0
for (( i=3; i<=$#; i++ ))
do
	arg=${!i}
	n=$((i+1))
	narg=${!n}
	if [ ${arg:0:1} == "-" ]; then
		if [ ${arg:1:1} == "-" ]; then
			if [ $arg == "--scratch" ]; then
				scratch=$narg
				delete_scratch=1
				lastok=1
			elif [ $arg == "--just-analyze" ]; then
				justanalyze=1
			elif [ $arg == "--find-parameters" ]; then
				findpars=1
			elif [ $arg == "--delete-failed" ]; then
				deletefailed=1
			elif [ $arg == "--runs" ]; then
				runs=$narg
				args="$args $arg $narg"
				lastok=1
			elif [ $arg == "--minutes" ]; then
				minutes=$narg
				args="$args $arg $narg"
				lastok=1
			elif [ $arg == "--quiet" ]; then
				quiet=1
				args="$args $arg"
			else
				args="$args $arg"
				lastok=1
			fi
		else
			arg=${arg:1}
			len=${#arg}
			if [ $len -eq 0 ]; then
				usage
			fi
			for (( j=0; j<$len; j++ ))
			do
				option=${arg:0:1}
				if [ $option == "J" ]; then
					justanalyze=1
				elif [ $option == "F" ]; then
					findpars=1
				elif [ $option == "D" ]; then
					deletefailed=1
				elif [ $option == "q" ]; then
					quiet=1
					args="$args -q"
				else
					lastok=1
					if [ $j -eq $((len-1)) ]; then
						if [ $option == "S" ]; then
							scratch=$narg
							delete_scratch=1
							lastok=2
						elif [ $option == "r" ]; then
							runs=$narg
							args="$args -r"
						elif [ $option == "m" ]; then
							minutes=$narg
							args="$args -m"
						else
							args="$args -$option"
						fi
					else
						usage
					fi
				fi
				arg=${arg:1}
			done
		fi
	else
		if [ $lastok -eq 0 ]; then
			usage
		else
			if [ $lastok -eq 1 ]; then
				args="$args $arg"
			fi
			lastok=0
		fi
	fi
done

# iterate through the parameter sets file line by line
par=0
sets=`cat $parsets`

for line in ${sets[@]}; do
	# separate the mutable parameter values from the rest of the line
	her1=${line%%,*}
	line=${line#*,}
	her7=${line%%,*}
	line=${line#*,}
	her13=${line%%,*}
	line=${line#*,}
	delta=${line%%,*}
	line=${line#*,}
	
	# run simulations for this parameter set
	run_sim $her1 $her7 $her13 $delta $line
	
	# increment the parameter set number
	((par++))
done

