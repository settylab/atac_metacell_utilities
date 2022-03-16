rule all_seqs:
    input: 
        peaks = config["results"] + "peaks.bed"
        script = config["scripts"] + "seq_gl.R"
    params:
        span = 150
        org = "hg38"
    output: 
        outfile = "all_seqs.fa"
    shell:
        "Rscript input.script input.peaks output.outfile params.span params.org"


rule fimo:
    input: 
    	meme = config["in_meme"]
        seqs = "all_seqs.fa"
    params:
        outputdir = "results/fimo_result"
    output:
        directory("results/fimo_result")
    shell:
        fimo  -oc params.outputdir input.meme input.seqs
        
