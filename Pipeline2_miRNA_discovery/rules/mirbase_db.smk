rule make_mirbase_db:
    input:
        fasta=config["mirbase"]
    output:
        "reference/miRBase_db.nsq"
    shell:
        """
        makeblastdb -in {input.fasta} -dbtype nucl -out reference/miRBase_db
        """
