rule write_MethylMatrix:
    input:
        expand('references/reference_methylomes/singleBP_bedGraph/{methylome}.singlebp.bedGraph', methylome=config['REFERENCE_METHYLOMES'])
    output:
        MethylMatrix = 'references/reference_methylomes/MethylMatrix/MethylMatrix.txt',
        MethylMatrix_header = 'references/reference_methylomes/MethylMatrix/MethylMatrix.header',
        MethylMatrix_body = 'references/reference_methylomes/MethylMatrix/MethylMatrix.body'
    message: "Creating MethylMatrix"
    run:
        files_to_merge = ' '.join(natsorted(input))
        f = open(METHYLOME_TABLE)
        headers = ' '.join(natsorted([x.split('\t')[0] for x in f.readlines()[1:]]))
        f.close()
        print(files_to_merge)
        print(headers)
        unionbed_command = f'bedtools unionbedg -i {files_to_merge} -header -names {headers} -filler - > {output.MethylMatrix}'
        os.system(unionbed_command)
        body_only_command = f'tail -n +2 {output.MethylMatrix} > {output.MethylMatrix_body}'
        os.system(body_only_command)
        header_only_command = f'head -n 1 {output.MethylMatrix} > {output.MethylMatrix_header}'
        os.system(header_only_command)

rule split_MethylMatrix_by_chr:
    input:
        MethylMatrix_header='references/reference_methylomes/MethylMatrix/MethylMatrix.header',
        MethylMatrix_body = 'references/reference_methylomes/MethylMatrix/MethylMatrix.body'
    output:
        MethylMatrix_chr = 'references/reference_methylomes/MethylMatrix/{chr}/MethylMatrix_{chr}'
    shell:
        """
        cp {input.MethylMatrix_header} {output.MethylMatrix_chr}
        bedextract {wildcards.chr} {input.MethylMatrix_body} >> {output.MethylMatrix_chr}
        """
