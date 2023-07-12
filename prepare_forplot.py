#infile="/cluster/work/pausch/naveen/GEMMA/BIMBAM/hd/sire/milk_production/CHR25/assoc.txt"
#outfile ="test.txt"

infile = snakemake.input.infile
outfile = snakemake.output.outfile


out = open (outfile, "w")
with open (infile) as inf:
    for lnum,line in enumerate (inf):
        spl=line.rstrip().split()
        if lnum==0:
            header=spl
        else:
            info=dict (zip (header, spl))
            mychr, mypos =info["rs"].split ("_")
            out.write (f"{mychr}\t{mypos}\t{info ['af'] }\t{ info['p_wald'] }\t{ info['p_lrt'] }\t{ info['p_score'] }\n")
            
        
out.close ()
