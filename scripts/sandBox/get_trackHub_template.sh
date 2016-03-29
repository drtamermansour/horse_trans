track_hub="$1"

cd $track_hub
wget -r --no-parent --reject "index.html*" -nH --cut-dirs=3 http://genome.ucsc.edu/goldenPath/help/examples/hubDirectory/
rm robots.txt
mv hubDirectory/* .
rm -R hubDirectory