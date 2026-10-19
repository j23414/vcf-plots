# vcf-plots

Exploratory visualization of VCF datasets, originally focused on influenza A virus samples but these scripts will be generalized as more projects arrive.

<img src="https://raw.githubusercontent.com/j23414/vcf-plots/refs/heads/main/diagram.png" width="350">

Expected inputs:

* VCF - variants in wide format, expecting exact IAV specific chromosomes
* Depth - samtools depth for all positions in a chromosome, to evaluate the SNPs

Exporatory visualizations in either static in a local Shiny Dashboard

* SNPs along the genome plots where y=minorPercentage
* Depth along the genome, normalize the minorPercentage based on the middle half of the depth plot

