#qsub -l mem_requested=2G -t 1-4 -S /bin/bash -cwd -o /net/eichler/vol8/home/nkrumm/EXOMES/zeromq_test/mapper/SGE_logs -e /net/eichler/vol8/home/nkrumm/EXOMES/zeromq_test/mapper/SGE_logs -j y mapperStart.sh 

source ~/.bash_profile
module load modules modules-init modules-gs modules-eichler

python mrsfast_wrapper.py e171
