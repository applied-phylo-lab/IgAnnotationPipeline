import glob
import yaml
import platform
import os
from pathlib import Path

OUTPUT="/local/storage/dhardesty/assemblies/#Wall_Lizards/"
mapping="LatinToCommonToFilename.csv"
IgDetective_dir="/local/storage/kav67/IgDetective/"
Human_ref="Human_ref/"

loci = ["IGH","IGL","IGK","TRA","TRB","TRG"] 

files = [
    f"{p.parent.name}/{p.stem}"
    for p in Path(OUTPUT).glob("*/*")
    if p.is_file() and (p.name.endswith("pri.fna") or p.name.endswith("alt.fna"))
]

print(files)

rule all:
    input:
        expand(OUTPUT + "{file}/{l}_table.tsv", file=files, l=loci)

#comment out the following rule if you want to skip the IgDetective step
rule IgDetective:
    resources:
        mem="2G"
    threads: 10
    log: os.path.join(OUTPUT,"logs", "IgDetective_{file}_{l}.log")
    input:
        fasta = OUTPUT + "{file}.fna"
    output:
        csv = OUTPUT + "{file}/combined_genes_{l}.txt"
    params:
         dir = OUTPUT + "{file}/",
         name = mapping,
         IgDetective = IgDetective_dir
    shell:
        "echo " + platform.node() + " >> {log} 2>&1 && \
        cd {params.IgDetective} &&  \
        python run_iterative_igdetective.py {input.fasta} {params.dir} >> {log} 2>&1"


rule createDBfromIgDetective:
    resources:
        mem="2G"
    threads: 10
    log: os.path.join(OUTPUT,"logs", "createDBfromIgDetective_{file}_{l}.log")
    input:
        script = 'csvToFasta.py',
        csv = OUTPUT+"{file}/combined_genes_{l}.txt",
        mapping = mapping
    output:
        fasta = OUTPUT + "{file}/{l}V.fasta"
    params:
         dir = OUTPUT + "{file}",
         locus = "{l}"
    shell:
        "echo " + platform.node() + " >> {log} 2>&1 && \
        python {input.script} {input.csv} {params.dir} {input.mapping} {params.locus}>> {log} 2>&1"

#uses human as default species
rule vquest:
    resources:
        mem="10G",
    threads: 10
    log: os.path.join(OUTPUT, "logs", "findVDJ_{l}_{file}.log")
    input:
        script = 'IMGT_vquest.py',
        fasta = OUTPUT + "{file}/{l}V.fasta"
    output:
        vquest = OUTPUT + "{file}/{l}_vquest/vquest_airr.tsv"
    params:
        dir = OUTPUT + "{file}/{l}_vquest/",
        locus = "{l}"
    shell:
        "echo " + platform.node() + " &>> {log} && \
         python {input.script} -f {input.fasta} -o {params.dir} -l {params.locus} -s human &>> {log}"

rule create_v_ref_gapped:
    resources:
        mem="2G"
    threads: 10
    log: os.path.join(OUTPUT, "logs", "create_v_ref_gapped_{l}_{file}.log")
    input:
        script = 'create_v_ref_gapped.py',
        vquest = OUTPUT + "{file}/{l}_vquest/vquest_airr.tsv"
    output:
        v_ref_gapped = OUTPUT + "{file}/{l}V_gapped.fasta"
    params:
        dir = OUTPUT + "{file}/",
        locus = "{l}"
    shell:
        "echo " + platform.node() + " &>> {log} && \
         python {input.script} -i {input.vquest} -o {params.dir} -l {params.locus} &>> {log}"


rule digger:
    resources:
        mem="10G"
    threads: 10
    log: os.path.join(OUTPUT, "logs", "digger_{l}_{file}.log")
    input:
        fasta = OUTPUT + "{file}.fna",
        vref = OUTPUT + "{file}/{l}V.fasta",
        vref_gapped = OUTPUT + "{file}/{l}V_gapped.fasta"
    output:
        combined = OUTPUT + "{file}/{l}_ref.fasta",
        digger = OUTPUT + "{file}/{l}_digger_r.tsv",
        diggerf = OUTPUT + "{file}/{l}_digger_f.tsv"
    params:
        locus = "{l}",
        dir = OUTPUT + "{file}/",
        jref = lambda wildcards: 
            os.path.join(OUTPUT, f"{wildcards.file}/{wildcards.l}J.fasta") 
            if os.path.exists(os.path.join(OUTPUT, f"{wildcards.file}/{wildcards.l}J.fasta")) 
            else os.path.join(Human_ref, f"Homo_sapiens_{wildcards.l}J.fasta"),
        dref = lambda wildcards: (
            os.path.join(OUTPUT, f"{wildcards.file}/{wildcards.l}D.fasta") 
            if wildcards.l in ["IGH", "TRB"] and os.path.exists(os.path.join(OUTPUT, f"{wildcards.file}/{wildcards.l}D.fasta"))
            else os.path.join(OUTPUT, "Human_ref", f"Homo_sapiens_{wildcards.l}D.fasta") 
            if wildcards.l in ["IGH", "TRB"] 
            else None
        )
    shell:
        """
        echo $(hostname) &>> {log} && \

        if [ -f {params.dref} ] && [ ! -z "{params.dref}" ]; then
            cat {input.vref} {params.dref} {params.jref} > {output.combined}
            DREF_OPTION="-d_ref {params.dref}"
        else
            cat {input.vref} {params.jref} > {output.combined}
            DREF_OPTION=""
        fi

        cd {params.dir}

        digger {input.fasta} \
        -v_ref {input.vref} \
        $DREF_OPTION \
        -j_ref {params.jref} \
        -v_ref_gapped {input.vref_gapped} \
        -ref igd,{output.combined} \
        -species human \
        -sense reverse \
        -locus {params.locus} \
        {output.digger} &>> {log}

        digger {input.fasta} \
        -v_ref {input.vref} \
        $DREF_OPTION \
        -j_ref {params.jref} \
        -v_ref_gapped {input.vref_gapped} \
        -ref igd,{output.combined} \
        -species human \
        -sense forward \
        -locus {params.locus} \
        {output.diggerf} &>> {log}
    """

rule combine_results:
    resources:
        mem="2G"
    threads: 10
    log: os.path.join(OUTPUT, "logs", "combine_results_{l}_{file}.log")
    input:
        script = 'combine_results.R',
        IgDetective = OUTPUT + "{file}/combined_genes_{l}.txt",
        digger = OUTPUT + "{file}/{l}_digger_r.tsv",
        diggerf = OUTPUT + "{file}/{l}_digger_f.tsv"
    output:
        combined = OUTPUT + "{file}/{l}_table.tsv",   
    params:
        dir = OUTPUT + "{file}/",
        locus = "{l}"
    shell:
        "echo " + platform.node() + " &>> {log} && \
         Rscript {input.script} --base_dir {params.dir} -l {params.locus} &>> {log}"
