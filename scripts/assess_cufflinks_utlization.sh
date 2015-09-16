sample_list="$1"
cufflinks_utlization="$2"

echo $(pwd) >> $cufflinks_utlization
while read dir; do
  f=$dir/cufflinks.e*
  if [ -f $f ]; then
    f2=$(echo $f | sed 's/cufflinks.e/cufflinks.o/')
    echo $dir >> $cufflinks_utlization
    grep "resources_used.vmem =" $f2 >> $cufflinks_utlization
    grep "resources_used.walltime =" $f2 >> $cufflinks_utlization
    grep "Resource_List.nodes =" $f2 >> $cufflinks_utlization
  else echo "no cufflinks report in $dir" >> $cufflinks_utlization
fi; done < $sample_list