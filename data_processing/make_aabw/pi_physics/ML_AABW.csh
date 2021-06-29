#!/bin/bash -l
#SBATCH -J ML_MOCsig_process
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10
#SBATCH -n 1
#SBATCH -t 23:00:00
#SBATCH -A NCGD0039 
#SBATCH -p dav
#SBATCH -e ML_MOCsig_process.err.%J
#SBATCH -o ML_MOCsig_proces_processs.out.%J

export TMPDIR=/glade/scratch/duvivier/moc-temp

module load ncl

export casename=b.e21.B1850.f09_g17.CMIP6-piControl.001b
#export casename=b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b
export outdir=/glade/scratch/duvivier/CESM2/MOC/
export indir=/glade/campaign/cgd/ppc/dbailey/${casename}/ocn/proc/tseries/month_1/
export timestamp="092101-097012"
scriptdir=/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/analysis/data_processing/make_aabw/
#scriptdir=`pwd`

for ((i=921; i<971; i++)); do 
 export year=$i
 echo "working on: " $i
 ncl ${scriptdir}/MOCeul.sigmacoord.POP.mymocloops.timeseries.ncl
done

echo "Finished!"


