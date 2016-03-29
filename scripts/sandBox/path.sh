myRoot=$"/mnt/ls15/scratch/users/mansourt/Tamer"
source $myRoot/config.txt

cufflinks_run="nonGuided_Cufflinks"
cuffmerge_run="nonGuided_Cuffmerge"
isoformfrac=0.05
tissue_Cuffmerge=$tissue_merge/cuffmerge
dist_dir="all_tissues_frac$isoformfrac"
cuffmerge_output=$dist_dir/$cufflinks_run/$cuffmerge_run

