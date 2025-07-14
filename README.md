# IgAnnotationPipeline
Pipeline for IG Annotation with IgDetective &amp; Digger

To install snakemake pipeline:
```
git clone https://github.com/applied-phylo-lab/IgAnnotationPipeline.git
cd IgAnnotationPipeline

conda env create -f annotation.yml
conda activate annotation
```

To run snakemake pipeline:
Change the following lines in the Snakefile:
```
OUTPUT="/path/to/results"
IgDetective_dir="/path/to/IgDetective_directory"
```

Also change the contents of LatinToCommonToFilename.csv according to your samples

Once everything is setup run snakemake:
```
# First test if there are any errors
snakemake -n

# If no errors show up run:
snakemake --jobs 100 #change according to availability of cores
```
