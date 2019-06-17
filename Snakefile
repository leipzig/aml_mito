import glob

FASTA = glob.glob('output/*.fa')
MMS = [fa+'.mm.py17.txt' for fa in FASTA]

rule default:
        input: "newseq_output.fa", MMS

rule catmms:
        input: MMS
        output: "allmms.txt"
        shell:
                """
                head -n 1 {input[0]} > {output}
                cat {input} | grep -v 'query' >> {output}
                """

rule reseq:
        input: "reseq.pl"
        output: "newseq_output.fa"
        shell:
                """
                rm -f mitoout.xls newseq_output.fa recalls_output.txt
                perl reseq.pl
                """

rule mitomaster:
        input: "{file}.fa"
        output: "{file}.fa.mm.py17.txt",
        shell:
                """
                perl mitomasterclient.pl {input} > {output}
                """

rule clean:
       shell:
                """
                rm -f output/*
                rm -f mitoout.xls newseq_output.fa recalls_output.txt
                """
