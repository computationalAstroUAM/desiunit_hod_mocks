#! /bin/bash
cd ${PBS_O_WORKDIR}

indir=${indir}
outdir=${outdir}
red=${red}
what2tar=${what2tar}

tar -zcvf ${outdir}${red}/${red}_${what2tar}.tar.gz ${indir}${red}/*_${what2tar}_*
