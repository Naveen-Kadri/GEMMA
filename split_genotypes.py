import os
infile=snakemake.input.infile
prefix =snakemake.params.prefix
nvar =snakemake.params.nvar

if not os.path.exists (prefix):
    os.makedirs (prefix)
    
part_number = 1
n=0
with open (infile) as inf:
    for i, line in enumerate(inf):
        n+=1
        if n == 1:
            out= open (prefix + "part_" +   str(part_number) + ".txt", "w")
        out.write (f"{line}")

        if n==nvar:
            out.close()
            n=0
            part_number+=1
