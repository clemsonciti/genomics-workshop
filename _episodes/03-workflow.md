---
title: "Variant calling workflow"
teaching: 20
exercises: 1
---

# Downloading Data
## Fetch from SRA database

~~~
prefetch -v SRR527257
~~~

## Convert SRA to fastq	

~~~
fastq-dump $file --outdir ~/data_carpentry/
~~~

