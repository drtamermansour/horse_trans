The Horse Transcriptome project 
===============================
Most of the gene models of the current horse genome lack trancriptional evidence support. This project develops a pipeline for analysis of RNAseq to fullfill these aims:

1. Make use of RNAseq as a transcriptional evidence to provide more accurate gene models.
2. Compare tissue specific gene models
3. Test different approaches for effective integrative analysis of several RNAseq experiments 
4. Use GitHub and UCSC genome browser to provide a scientific collaborative platform for laboratories interested in the horse functional genomics and provide a replicable model for other organisms. 
5. Allow continuous update of transcriptomes with new RNAseq experiments. Also facilitate reproduction of the whole annotation pipeline with newer versions of the genome 
   
This github project contains all the scripts required to reproduce the analysis of the horse transcriptome project on the HPC of Michigan state university. Also it has the basic tree structure of the UCSC track hub. This design will allow the community to share in developiong the track hub and also they can fork to have their own version with different view options.

Coding guidelines in this project
---------------------------------
- Fork the horse_trans repository to your Github account, clone it to your local machine, and Configure Git to sync your fork with the original repository (if you want). To learn more about doing this, you can check the [github documenation](https://help.github.com/articles/fork-a-repo/)  
- The "main.sh" script represents the index page of the whole project. The first section of the main script should construct the basic diretory structure for you.
- Preparation of the input data:
      <p align="center">
         <img src="directory_structure.png" width="300"/>
      </p>
   * Every tissue has a separate folder carrying its name (maximum 14 letter). All tissue folders should be housed in the prepdata folder (a folder that was created in the previous step)
   * Every RNAseq library should have a separate folder in the corresponding tissue folder.
   * The library folder name should have this format:
      - start with PE_ or SE_ according to the sequencing type
      - Then it should should have the read length followed by underscore
      - Then it should have fr.unstranded_ , fr.firststrand_ , fr.secondstrand_ according to lib type
      - Then the owner name (or names separted by dots) followed by underscore
      - Then the date of sequencing as MMDDYYYY
   * The raw data files should be kept in a folder named fastq_data in the library folder
   * The raw data files should fullfill these criteria 
      - The file names should fit the format \*\_R1_\*.fastq.gz & \*\_R2_\*.fastq.gz for PE reads or \*\_SR_\*.fastq.gz for SE
      - All reads in every given file should belong to ONE sequencing lane.
      - If there are sample replicates from different lanes, you can add a text file called "replicates.txt". Each line in this file should have the names of one sample replicates with space separation. (only the \*\_R1_\*.fastq.gz for PE reads or \*\_SR_\*.fastq.gz for SE)
      - The first syllabus (the part of the file name before \_R1_ , \_R2_ or \_SR_) should be unique
      - All samples should be prepared so that they have enconding "Sanger / illumina 1.9"
   * prep_proj.sh script can automatically convert an SRA repository into the appropriate format (the script assumes no sample replicates)
- All the analysis was done on the High performance computer of Michigan state university. The main script represents an abstraction administrator. It does not run any of the programs but always initiates daughter scripts to do the job. This allows new forks or branches to provide different implementations of the daughter scripts to run the same pipeline with different versions or on different machines.
- With such a large scale analysis, it is difficult to check that every job ended successfully. I tried to add a check point after every step to make sure that I am getting out what I really expect. I believe adding more of these tests would be a good practice
- Every step in the analysis receives the input data as a sample list. This enables the users to run the pipeline for selected samples. Also enables easier re-running of failed samples. On the long term this would allow the community to keep the transcriptome list updated without re-analysis or re-coding. 
- Track hubs: 
   * Each hub should include tracks for set of samples that underwent the same "non-branching" analytical pipeline. 
   * Each hub should contain tracks for every tissue. Tissues with mutliple libraries would have a special composite track which allow the visualization of the separate libraries. The hub also should contain track that represent the integration of all tissues.  
   * I tried to automate the creation and edition of the tracks so that the user just pass a list of all the annotation files or even only the newer files to the "edit_trackDb.sh" script which takes care of arranging the libraries and tissues into the appropriate format. Hopefully it is going to work the way I like.

Current avaliable track hubs
============================
You can add them to your hubs on UCSC using these URLs:

-  Tophat2 followed by Non guided Cufflinks/Cuffmerge:
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopNonGuidCuff.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopNonGuidCuff.txt)
-  Tophat2 followed by Non guided Cufflinks/Cuffmerge (The unfiltered version):
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TNGCuffUnfilt.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TNGCuffUnfilt.txt) 
-  Adding UTRs to the coding transcripts:
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopCuffmerge_TD.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_TopCuffmerge_TD.txt) 
-  Publically avaliable assemblies
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_public_assemblies.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_public_assemblies.txt) 

<!---  
-  Diginorm followed by refGTFguided Tophat2 then refGTFguided Cufflinks https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_rawdigi_TopCuff.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_rawdigi_TopCuff.txt)
[//]: #    -  Assemblies after exom merge
https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_exonMerge_assemblies.txt or you can view it directly by following this [link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=equCab2&hubUrl=https://raw.githubusercontent.com/drtamermansour/horse_trans/master/track_hub/hub_exonMerge_assemblies.txt)
--->
To add a hub, go to the genome browser home bage. In the My Data menu, open the Track Hubs page and click the [My Hubs](https://genome.ucsc.edu/cgi-bin/hgHubConnect) tab. This tab lists the unlisted track hubs that you have loaded into your browser. To import a new hub, type its URL into the text box, then click the Add Hub button.

Links for Downloads:
--------------------
-  [GTF file of filtered assembly](http://de.iplantcollaborative.org/dl/d/DFB23E2B-F749-4AC9-99BB-BC05BEDFAF79/filtered_Alltissues_Assembly.GTF)
-  [Transcriptome edited by common variants](http://de.iplantcollaborative.org/dl/d/5F00CC13-5775-4AD5-A814-E3AFC768E2D9/varFixed_Transcriptome.fa)
-  [Tabulated gene TPM expression for all tissues](http://de.iplantcollaborative.org/dl/d/8AD5668A-02C3-4E54-AEE3-509E15C54594/allTissues_geneTPM)  
-  [Tabulated isoform TPM expression for all tissues](http://de.iplantcollaborative.org/dl/d/FA197031-71D3-4AAD-911A-7ACDF7516911/allTissues_isoformTPM)  

