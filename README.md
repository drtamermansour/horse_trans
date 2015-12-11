The Horse Transcriptome project 
===============================
Most of the gene models of the current horse genome lack trancriptional evidence support. This project develops a pipeline for analysis of RNAseq to fullfill these aims:

1. Make use of RNAseq as a trascriptional evidence to provide more accurate gene models.
2. Compare tissue specific gene models
3. Test different approaches for effective intgrative analysis of several RNAseq experiments 
4. Use GitHub and UCSC genome browser to provide a scientific collaprative plateform for laboratories intrested in the horse functional genomics and provide a replicable model for other organisms. 
5. Allow continuous update of transcriptomes with new RNAseq expermints. Also facilitate reproduction of the whole annoatation pipeline with newer versions of the genome 
   
This github project contains all the scripts required to reproduce the analysis of the horse transcriptome project on the HPC of Michigan state university. Also it has the basic tree structure of the UCSC track hub. This design will allow the community to share in developiong the track hub and also they can fork to have their own version with different view options.

Coding guidelines in this project
---------------------------------
- All the scripts of the project are housed in the script folder
- The "main.sh" script represents the index page of the whole project. Scripts that I used for any special functions and not referenced in the main script are found in a subfolder called "special".
- All the directory structure of the projected goes back to a common working directory. You need to assess the path of your working directory to the variable myRoot at the beginning of the main script
- Preparation of the input data:
   * Every tissue has a separate folder carrying its name (maximum 14 letter)
   * Every RNAseq libarary should have a separate folder in the corresponding tissue folder.
   * The libarary folder name should have this format:
      - start with PE_ or SE_ according to the sequencoing type
      - Then it should should have the read length followed by underscore
      - Then it should have fr.unstranded_ , fr.firststrand_ , fr.secondstrand_ according to lib type
      - Then the owner name (or names separted by dots) followed by underscore
      - Then the date of sequencing as MMDDYYYY
   * The raw data files should be kept in a folder named fastq_data in the libarary folder
   * The raw data files should fullfill these criteria 
      - The file names should fit the format \*\_R1_\*.fastq.gz & \*\_R2_\*.fastq.gz for PE reads or \*\_SR_\*.fastq.gz for SE
      - All reads in every given file should belong to ONE sequencing lane.
      - If there are sample replicates from different lanes, you can add a text file called "replicates.txt". Each line in this file should have the names of one sample replicates with space separation. (only the \*\_R1_\*.fastq.gz for PE reads or \*\_SR_\*.fastq.gz for SE)
      - The first syllbus (the part of the file name before \_R1_ , \_R2_ or \_SR_) should be unique
      - All samples should be prepared so that they have enconding "Sanger / illumina 1.9"
   * prep_proj.sh script can autmatically convert an SRA reposatory into the appropriate format (the script assume no sample replicates)
- All the analysis was done on the High performance computer of Michigan state university. The main script represents an abstarction adminstroter. It does not run any of programs but always initiate daughter scripts to do the job. This allows new forks or branches to provide different implementations of the daughter scripts to run the same pipeline with different versions or on different machines.
- With such a large scale analysis, it is difficult to check that every job ended successfully. I tried to add a check point after every step to make sure that I am getting out what I really expect. I believe adding more of these tests would be a good practice
- Every step in the analysis recive the input data as a sample list. This enable the users to run the pipeline for selected samples. Also enable easier re-runing of failed samples. On the long trem this would allow the community to keep updated the tarnscriptome list without re-analysis or re-coding. 
- Track hubs: 
   * Each hub should include tracks for set of samples that underwent the same "non-branching" analytical piepline. 
   * Each hub should contain tracks for every tissue. Tissues with mutliple libraries would have a special composite track which allow the visualization of the separate libraries. The hub also should contain track that represent the integration of all tissues.  
   * I tried to automate the creation and edition of the tracks so that the user just pass a list of all the annotation files or even only the newer files to the "edit_trackDb.sh" script which take care of of arranging the libraries and tissues into the appropriate format. Hopefully it is going to work the way I like.

Current avaliable track hubs
============================
You can add them to your hubs on UCSC using these URLs:

-  refGTFguided Tophat2 followed by Non guided Cufflinks/Cuffmerge:
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopNonGuidCuff.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopNonGuidCuff.txt)
-  refGTFguided Tophat2 followed by refGTFguided Cufflinks then Cuffmerge:
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopGuidedCuff.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopGuidedCuff.txt) 
-  Adding UTRs to the coding transcripts:
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopCuffmerge_TD.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopCuffmerge_TD.txt) 
-  Publically avaliable assemblies
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_public_assemblies.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_public_assemblies.txt) 
-  Diginorm followed by refGTFguided Tophat2 then refGTFguided Cufflinks
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_rawdigi_TopCuff.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_rawdigi_TopCuff.txt)
-  Assemblies after exom merge
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_exonMerge_assemblies.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_exonMerge_assemblies.txt)

To add a hub, go to the genome browser home bage. In the My Data menu, open the Track Hubs page and click the [My Hubs](https://genome.ucsc.edu/cgi-bin/hgHubConnect) tab. This tab lists the unlisted track hubs that you have loaded into your browser. To import a new hub, type its URL into the text box, then click the Add Hub button. 
