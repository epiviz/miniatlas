### Instructions for setting up Trackhubs

There are several ways to setup trachub instances. UCSC has instructions on setting up and creating trackhubs - <http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html#Setup>

- Brian Herb uses the R script (`BICCN_MiniAtlas_Trackhub.R`) in this folder to setup trackhubs
    - requires `bedtools` or `ucsctools` installed
- Use the `trackhub` python package - <https://github.com/daler/trackhub>


Always validate the hub instace by using the `hubCheck` utility from UCSC (<https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html#Compatibility>)