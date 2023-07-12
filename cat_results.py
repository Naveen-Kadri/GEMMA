infiles=snakemake.input.infiles
out= open (snakemake.output.outfile, "w")

print (f"infiles are {infiles}")


for fnum,myfile in enumerate(infiles):
    with open (myfile,"rt") as inf:
        for lnum,line in enumerate(inf):
            if fnum>0 and lnum==0:
                continue
            out.write (f"{line}")
out.close()
            
