rule mirdeep2:
    input:
        reads=f"{PROJECT}/filtered/{{sample}}.collapsed.fa",
        arf=f"{PROJECT}/filtered/{{sample}}.arf",
        genome=config["genome_fasta"]
    output:
        done=f"{PROJECT}/mirdeep/{{sample}}/result_done.txt",
        report=f"{PROJECT}/mirdeep/{{sample}}_report.txt"
    conda:
        "../envs/mirdeepP2.yaml"
    shell:
        """
        OUTDIR={PROJECT}/mirdeep/{wildcards.sample}
        RUN_DIR=$OUTDIR/run
        mkdir -p $RUN_DIR
        cd $RUN_DIR
        miRDeep2.pl \
            {input.reads} \
            {input.genome} \
            {input.arf} \
            none none none \
            -t plant -d
        cp result_*.csv {output.report}
        touch {output.done}
        """
