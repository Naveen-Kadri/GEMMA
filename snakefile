from collections import defaultdict
tfile='trait_info.txt'

tnumbers = defaultdict (dict)
cats =defaultdict (list)
with open (tfile) as inf:
    for lnum,line in enumerate (inf):
        if lnum==0:
            continue
        trait, ocat, tdef, cat=line.rstrip().split ("\t")
        cat=cat.replace (" ", "_")
        cats [cat].append (trait)


for sex in ["sire", "dam"]:
    tnumber = dict ()
    phenofile =f"/cluster/work/pausch/naveen/GEMMA/BIMBAM/hd/{sex}/CHR25/traits.txt"
    with open (phenofile) as inf:
        for tnum,line in enumerate(inf):
            trait=line.rstrip()
            tnumber [trait] = tnum +1

    for cat in cats:
        mycat = cats [cat]
        mynumbers = []
        for trait in mycat:
            if trait in tnumber:
                mynumbers.append ( str(tnumber [trait])   )
        mynumbers=" ".join (mynumbers)
        tnumbers [sex] [cat] = mynumbers
    tnumbers[sex] ["all"] = " ".join (str(num)  for num in range(1,tnum+2))

#retrive columbers as tnumbers [sex] [cat] where cat are
print (tnumbers)




GEMMA ="/cluster/home/nkadri/PROGRAMS/GEMMA/gemma-0.98.5-linux-static-AMD64 "

##VARIABLES/WILDCARDS/FILES
chromosomes =range (1,30)
densities = ['hd', 'seq']
sexes = ['sire', 'dam']
##categories = ["Fertility", "milk_production", "Linear_conformation_trait", "all"]
#categories = ["Fertility", "milk_production", 'Teat']
categories = ['milk_production', 'Teat']
OUT_DIR = "/cluster/work/pausch/naveen/GEMMA/"
vcf_files={
    "seq" : "/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/combined.vcf.gz",
    "hd"  : "/cluster/work/pausch/naveen/SNPDATA/hd_imputed/CHR{chr}/combined.vcf.gz"
}

##prepared from files used earlier for GWAS [extreme values removed]
phenotype_file=OUT_DIR + "PHENOTYPES/{sex}/phenotypes.txt" 


rule all:
    input:
        expand(OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/assoc_new.txt", chr=chromosomes, sex=['sire', 'dam'], cat=['Teat', 'milk_production', 'Linear_conformation_trait'], density=['seq'])
        

#phenotype file is written for every chromosome -- it is ok for now
rule get_bim_bam:
    input:
        vcf_file=lambda wildcards : vcf_files [wildcards.density],
        phenotype_file=phenotype_file
    output:
        genotypes=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/genotypes.txt",
        phenotypes=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/phenotypes.txt",
        ids=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/ids.txt",
        traits=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/traits.txt"
    params:
        r2_thresh=0.3
    resources:
        mem_mb=16000,
        walltime="10:00"
    script:
        "get_bim_bam.py"


checkpoint split_genotypes:
    input:
        infile=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/genotypes.txt"
    output:
        directory(OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/parts/"),
        touch (OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/parts/splitting.done")
    resources:
        mem_mb=4000,
        walltime="04:00"
    params:
        prefix=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/parts/",
        nvar=50_000
    script:
        "split_genotypes.py"    
        
##MAYBE JUST USE HD snps?? [rules to extract the HD snps and combine into one file]

rule combine:
    input:
        infiles=expand (OUT_DIR + "BIMBAM/{{density}}/{{sex}}/CHR{chr}/genotypes.txt",chr=chromosomes)
    output:
        outfile=OUT_DIR + "BIMBAM/{density}/{sex}/GRM/genome_genotypes.txt"
    resources:
        mem_mb=4000,
        walltime="04:00"
    shell:
        "cat  {input.infiles} > {output.outfile}"    


#type 1 is for centered, authors say this provides better control for pop strat.
rule get_grm:
    input:
        genotypes = OUT_DIR + "BIMBAM/{density}/{sex}/GRM/genome_genotypes.txt",
        phenotypes=OUT_DIR + "BIMBAM/{density}/{sex}/CHR25/phenotypes.txt",
    output:
        log=OUT_DIR + "BIMBAM/{density}/{sex}/GRM/grm.log.txt",
        grm=OUT_DIR + "BIMBAM/{density}/{sex}/GRM/grm.cXX.txt"
    params:
        wd=OUT_DIR + "BIMBAM/{density}/{sex}/GRM/",
        prefix="grm",
        type=1,
        output="z"
    threads:
        20
    resources:
        mem_mb=24000,
        walltime="04:00"
    shell:
        "cd {params.wd} \n" +
        GEMMA + " -g {input.genotypes} -p {input.phenotypes} -gk {params.type}  -o {params.prefix} \n" +
        " mv output/grm.cXX.txt {output.grm} \n" + 
        " mv output/grm.log.txt {output.log}"


        
rule multivar:
    input:
        genotypes=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/parts/part_{part}.txt",
        phenotypes=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/phenotypes.txt",
        grm=OUT_DIR + "BIMBAM/hd/{sex}/GRM/grm.cXX.txt",
        flag=OUT_DIR + "BIMBAM/{density}/{sex}/CHR{chr}/parts/splitting.done"
    output:
        res=OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/assoc_{part}.txt",
        log=OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/assoc_{part}_log.txt",
    params:
        wd = OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/",
        trait_numbers=lambda wildcards : tnumbers [wildcards.sex] [wildcards.cat],
        prefix="part_{wildcards.part}"
    threads:
        1
    resources:
        mem_mb=8000,
        walltime="24:00"
    shell:
        "cd {params.wd} \n" + 
        GEMMA + " -lmm 4 -g {input.genotypes} -p {input.phenotypes} -n {params.trait_numbers} -k {input.grm} -o {params.prefix} \n" +
        " mv output/{params.prefix}.assoc.txt {output.res} \n" +
        " mv output/{params.prefix}.log.txt {output.log}\n"


def list_files (wildcards):
    mydir=checkpoints.split_genotypes.get (cat=wildcards.cat,chr=wildcards.chr,sex=wildcards.sex, density=wildcards.density).output [0]
    parts=glob_wildcards (mydir + "/part_{part}.txt").part
    parts = sorted ([int(el) for el in parts])
    outfiles = expand(OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/assoc_{part}.txt", density=wildcards.density, sex=wildcards.sex, cat=wildcards.cat, chr=wildcards.chr, part=parts)
    return (outfiles)
    
rule merge_multivar:
    input:
        infiles=list_files
    output:
        outfile=OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/assoc_new.txt"
    script:
        "cat_results.py"

        
rule prepare_forplot:
    input:
        infile=OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/assoc.txt",
    output:
        outfile=OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/CHR{chr}/forplot.txt",
    resources:
        mem_mb=2000,
        walltime="04:00"
    script:
        "prepare_forplot.py"

rule manhattan:
    input:
        infiles=expand (OUT_DIR + "BIMBAM/{{density}}/{{sex}}/{{cat}}/CHR{chr}/forplot.txt", chr=chromosomes)
    params:
        maf_thresh=0.001,
        chromosomes=chromosomes
    resources:
        mem_mb=8000,
        walltime="04:00"
    output:
        plot_file=OUT_DIR + "BIMBAM/{density}/{sex}/{cat}/manhattan_{density}_{cat}_{sex}.tiff",
    script:
        "manhattan.R"
