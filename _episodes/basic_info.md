## Variant Calling
For the genomics part of this workshop, we will be running a **variant
calling pipeline**, which starts with several raw sequencing data
files from several samples, and ends with a single data file in VCF
format that contains the Single Nucleotide Polymorphism (SNP) data.

SNPs are any location in the genome where at least one of the samples has a
mutation that gives it a different nucleotide than the other samples.

SNPs are used in many downstream applications in biology, including:
+  Genome-Wide Assocation Studies
+  Disease Risk Assessment
+  Tests of Relatedness between individuals
+  Population Genetics (tests of evolutionary questions)
+  Phylogenetics (constructing species trees)
+  Marker-assisted breeding

In today's exercise, after we generate our VCF file, we will use it to
plot a small phylogenetic tree, showing the relationship between our samples.

