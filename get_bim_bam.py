import gzip
# mychr=12
# vcf_file=f"/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{mychr}/combined.vcf.gz"
# phenotype_file="/cluster/work/pausch/naveen/GEMMA/PHENOTYPES/sire/phenotypes.txt"
# out_genotypes=open ("genotypes.txt", "w")
# out_phenotypes=open ("phenotypes.txt", "w")
# out_traits=open ("traits.txt", "w")
# out_ids=open ("ids.txt", "w")
# r2_thresh=0.3

vcf_file=snakemake.input.vcf_file
phenotype_file=snakemake.input.phenotype_file
out_phenotypes = open (snakemake.output.phenotypes,"w")
out_genotypes = open (snakemake.output.genotypes,"w")
out_ids=open (snakemake.output.ids, "w")
out_traits=open(snakemake.output.traits, "w")
r2_thresh=snakemake.params.r2_thresh



#hash phenotyped ids
phenotyped=dict()
with open (phenotype_file) as inf:
    for lnum,line in enumerate(inf):
        spl=line.rstrip().split("\t")
        if lnum==0:
            for el in spl[1:]:
                out_traits.write (f"{el}\n")
            out_traits.close()
        phenotyped [  spl[0] ] = "\t".join (spl [1:])
print(f"Number of phenotyped sires = {len  (phenotyped)  }")

        

genotype_order = []

all_var = 0
kept_var = 0
with gzip.open (vcf_file,"rt") as inf:
    for line in inf:
        if line[0:2]!="##":
            spl=line.rstrip().split("\t")
            if line[0:6] =="#CHROM":
                ids=spl[9:]
            else:
                all_var +=1
                #if all_var >1000:
                    #break
                if "," in spl[3] or "," in spl[4]:
                    continue
                infos = spl [7].split(";")
                if len (infos) ==3:
                    infos=infos [1:]
                try:
                    acc=float (infos[0].split("=") [1])
                except IndexError:
                    print (f"Index error at position {spl[1] }")
                    continue
                af=infos[1].split("=") [1]
                #print (infos,maf, acc)
                if acc < r2_thresh or af == "0" or af == "1":
                    continue
                kept_var +=1
                snpid = spl[0] + "_" + spl[1]
                gts=spl[9:]
                doses = []
                for myid,gt in zip(ids, gts):
                    if myid in phenotyped:
                        if kept_var == 1:
                            genotype_order.append (myid)
                        gt, dose=gt.split (":")
                        if dose=='.':
                            gdose=str (int (gt [0]) + int (gt [2]))
                            doses.append (gdose)
                        else:
                            doses.append(dose)
                tw=",".join([snpid, spl[3], spl[4]] + doses )
                out_genotypes.write (f"{tw}\n")
out_genotypes.close ()

print (f"Number of variants, Number kept  == {all_var}, {kept_var}")

for myid in genotype_order:
    out_phenotypes.write (f"{phenotyped[myid]}\n")
out_phenotypes.close()
out_ids.close()




