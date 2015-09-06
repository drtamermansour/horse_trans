UCSCgenome="$1"
hub_name="$2"
shortlabel="$3"
longlabel="$4"
track_hub="$5"


## create directory structure
mkdir -p $track_hub/$UCSCgenome/BigBed
> $track_hub/hub.txt
> $track_hub/genomes.txt
> $track_hub/$UCSCgenome/trackDb.txt
rm -f $track_hub/current_tracks

## create the hub file
echo "hub $hub_name" >> $track_hub/hub.txt
echo "shortLabel $shortlabel" >> $track_hub/hub.txt
echo "longLabel $longlabel" >> $track_hub/hub.txt
echo "genomesFile genomes.txt" >> $track_hub/hub.txt
echo "email drtamermansour@gmail.com" >> $track_hub/hub.txt


## create the genomes file
echo "genome $UCSCgenome" >> $track_hub/genomes.txt
echo "trackDb $UCSCgenome/trackDb.txt" >> $track_hub/genomes.txt




