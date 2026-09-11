"""
Metabarcoding pipeline based on EMU.
To run:
on thor: source /scratch/thdelcourt/metabarcoding/20250704/scripts/virtenv/EMU_ppl_py3.10/bin/activate
snakemake --cores 5 -s ../../../scripts/TDscripts/projects/metabarcoding/EMU_pipeline/pipeline.smk --configfile ./config_template.yaml --config workdir=/scratch/thdelcourt/metabarcoding/20250704/20250705_test_emu/test_ppl/Hoodia_gordonii_line1_3 sample_name=Hoodia_gordonii_line_1_20240813_1604_X1_ARW146_b6037b07 minQ=15 --use-conda
requires a virtual environment containing:
- snakemake
- SeqIO
- tabulate
- plotly
- numpy
"""



# ITS2_barcode_positions = {
#     79: "T",
#     87: "A",
#     122: "C",
#     126: "G",
#     158: "A"}

# dict_species_to_group = {'41737': 1, '41738': 1, '22217': 2, '22216': 2, '22218': 3}
# dict_groups_values = {1: ['41737', '41738'], 2:['22217', '22216'], 3:['22218']}
# dict_species_to_group = {'41738': 1}
# dict_groups_values = {1: ['41738']}

# Load config variables
workdir : config["workdir"]

sample_name = config["sample_name"]
minQ = config["minQ"]
print(minQ)
in_fastq_path = f"{config['in_fastq_dir']}/{sample_name}.fastq"
EMU_DB_path = config["EMU_DB_path"]
Hoodia_ref = config["Hoodia_ref"]

min_abundance = config["min_abundance"]
min_median_rap = config["min_median_rap"]
min_real_depth = config["min_real_depth"]
max_diff_cons_ref = config["max_diff_cons_ref"]

path_to_ViralMSA = config["path_to_ViralMSA"]
path_to_Hoodia_SNP_barcode_ref = config["path_to_Hoodia_SNP_barcode_ref"]
primers_path = config["primers_path"]
ref_seqs_DB_path = config["ref_seqs_DB_path"]
sequences_dict = config["sequences_dict"]
taxonomy_pickle = config["taxonomy_pickle"]
median_len_ITS2 = config["median_len_ITS2"]
window_size = config["window_size"]
trimmed_length = config["trimmed_length"]
path_to_identical_species_DB = config["path_to_identical_species_DB"]
path_to_ITS2_barcode_positions = config["path_to_ITS2_barcode_positions"]

def get_ITS2_length_range(median_len_ITS2, window=50, trimmed_length=50):
    """calculate the ITS2 length range: 450bp +-50bp
    return: tuple(min len, max len)"""
    median_len_ITS2_post_trimming = median_len_ITS2 - trimmed_length
    min_len_ITS2 = median_len_ITS2_post_trimming - int(window/2)
    max_len_ITS2 = median_len_ITS2_post_trimming + int(window/2)
    return  min_len_ITS2, max_len_ITS2

min_len_ITS2, max_len_ITS2 = get_ITS2_length_range(median_len_ITS2=median_len_ITS2, window=window_size, trimmed_length=50) # BUG!!! trimmed length cannot be hardcoded! (although it has always been the same until now)

rule pre_filtQ_QC_nanoplot:
    input:
        in_fastq_path
    output:
        pickle="pre_filtQ_QC/NanoPlot-data.pickle"
    shell:
        "ml load nanoplot ; "
        "NanoPlot -o pre_filtQ_QC --fastq {input} --tsv_stats -t 10 --raw --store  --no_static"
#
# rule pre_filtQ_QC_seqkit:
#     input:
#         in_fastq_path
#     output:
#         txt="pre_filtQ_QC/seqkit_stats.txt"
#     threads: 1
#     shell:
#         "ml load seqkit ; "
#         "seqkit stats {input} > {output.txt} "

rule plot_pre_filtQ_QC:
    input:
        pickle_data_path = rules.pre_filtQ_QC_nanoplot.output.pickle
    output:
        html_Q_hist="pre_filtQ_QC/Q_hist.html",
        html_Q_cum_hist="pre_filtQ_QC/Q_cum_hist.html",
        html_L_hist="pre_filtQ_QC/L_hist.html",
        html_L_hist_restrict="pre_filtQ_QC/L_hist_zoom.html",
        read_counts_table = "pre_filtQ_QC/read_counts_table.pickle"

    run:
        import plotly.graph_objects as go
        from plotly import plot
        from numpy import histogram, cumsum
        from pickle import load, dump
        import pandas

        # Load data
        with open(input.pickle_data_path, "rb") as in_stream:
            in_data = load(in_stream)

        # TODO: make more efficient by not creating dataframes. just count the booleans. not a major issue (from few tens of ms to a few ms)
        # TODO: add hardocded 50 (trimmed length) to trimmed_length defined in config
        # calculate reads counts and export to pickle
        # lengths + 50 to account for trimming downstream
        headers = ["", "All Q", f"Q > {minQ}"]
        out_list = []
        out_list.append((
        "reads in Total", len(in_data['lengths']), len(in_data['lengths'][(in_data['quals'] > minQ)])))
        out_list.append(("reads in range ITS2", len(
            in_data['lengths'][(in_data['lengths'] >= min_len_ITS2+50) & (in_data['lengths'] <= max_len_ITS2+50)]), len(
            in_data['lengths'][(in_data['lengths'] >= min_len_ITS2+50) & (in_data['lengths'] <= max_len_ITS2+50) & (
                        in_data['quals'] > minQ)])))

        with open(output.read_counts_table, "wb") as out_stream:
            dump(out_list, out_stream)

        # histogram of qualities
        series = in_data["quals"]
        nbins = 4 * round(max(series) - min(series))
        hist, bins = histogram(series, bins=nbins)
        fig = go.Figure(go.Bar(x=bins, y=hist,))
        fig.update_layout(xaxis_title="Qualities", yaxis_title="Count", bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.html_Q_hist, full_html=False, include_plotlyjs="cdn")

        # cumulative histogram of qualities
        series = in_data["quals"]
        hist, bins = histogram(series,bins=180)
        cum_hist = cumsum(hist)
        fig = go.Figure(go.Bar(x=bins,y=cum_hist,))
        fig.update_layout(xaxis_title="Qualities",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.html_Q_cum_hist, full_html=False, include_plotlyjs="cdn")

        # histogram of lengths
        bin_size = 100
        series = in_data["lengths"]
        nbins = round((max(series) - min(series)) / bin_size)
        hist, bins = histogram(series,bins=nbins)
        fig = go.Figure(go.Bar(x=bins,y=hist,))
        fig.update_layout(xaxis_title="Read length",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.html_L_hist, full_html=False, include_plotlyjs="cdn")

        # histogram of lengths with filtered out very long reads to facilitate histogram building at precision of 1bp
        series = in_data["lengths"].loc[lambda x: x <= 2500]
        nbins = max(series) - min(series)
        hist, bins = histogram(series,bins=nbins)
        fig = go.Figure(go.Bar(x=bins,y=hist,))
        fig.update_layout(xaxis_title="Read length",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.html_L_hist_restrict, full_html=False, include_plotlyjs="cdn")



rule filtQ:
    input:
        in_fastq_path
    output:
        "filtQ/out.fastq"
    shell:
        "ml load nanofilt ; "
        "NanoFilt -q {minQ} {input} > {output}"

rule post_filtQ_QC_nanoplot:
    input:
        rules.filtQ.output
    output:
        pickle="filtQ/nanoplot/NanoPlot-data.pickle"
    threads: 10
    shell:
        "ml load nanoplot ; "
        "NanoPlot -o filtQ/nanoplot/ --fastq {input} --tsv_stats -t 10 --raw --store  --no_static"

# rule post_filtQ_QC_seqkit:
#     input:
#         rules.filtQ.output
#     output:
#         txt="filtQ/seqkit_stats.txt"
#     threads: 1
#     shell:
#         "ml load seqkit ; "
#         "seqkit stats {input} > {output.txt} "

rule plot_post_filtQ_QC:
    input:
        pickle_data_path = rules.post_filtQ_QC_nanoplot.output.pickle
    output:
        html_Q_hist="post_filtQ_QC/Q_hist.html",
        html_Q_cum_hist="post_filtQ_QC/Q_cum_hist.html",
        html_L_hist="post_filtQ_QC/L_hist.html",
        html_L_hist_restrict="post_filtQ_QC/L_hist_zoom.html",
        read_counts_table = "post_filtQ_QC/read_counts_table.pickle"
    run:
        import plotly.graph_objects as go
        from numpy import histogram, cumsum
        from pickle import load, dump
        import pandas

        with open(input.pickle_data_path,"rb") as in_stream:
            in_data = load(in_stream)

        # calculate reads counts and export to pickle
        # lengths + 50 to account for trimming downstream
        headers = ["", "All Q", f"Q > {minQ}"]
        out_list = []
        out_list.append((
            "reads in Total", len(in_data['lengths']), len(in_data['lengths'][(in_data['quals'] > minQ)])))
        out_list.append(("reads in range ITS2", len(
            in_data['lengths'][
                (in_data['lengths'] >= min_len_ITS2 + 50) & (in_data['lengths'] <= max_len_ITS2 + 50)]), len(
            in_data['lengths'][
                (in_data['lengths'] >= min_len_ITS2 + 50) & (in_data['lengths'] <= max_len_ITS2 + 50) & (
                        in_data['quals'] > minQ)])))

        with open(output.read_counts_table,"wb") as out_stream:
            dump(out_list,out_stream)

        series = in_data["quals"]
        # histogram of qualities
        nbins = 4 * round(max(series) - min(series))
        hist, bins = histogram(series,bins=nbins)
        fig = go.Figure(go.Bar(x=bins,y=hist,))
        fig.update_layout(xaxis_title="Qualities",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.html_Q_hist, full_html=False, include_plotlyjs="cdn")

        # cumulative histogram of qualities
        hist, bins = histogram(series,bins=180)
        cum_hist = cumsum(hist)
        fig = go.Figure(go.Bar(x=bins,y=cum_hist,))
        fig.update_layout(xaxis_title="Qualities",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.html_Q_cum_hist, full_html=False, include_plotlyjs="cdn")

        # histogram of lengths
        bin_size = 100
        series = in_data["lengths"]
        nbins = round((max(series) - min(series)) / bin_size)
        hist, bins = histogram(series,bins=nbins)
        fig = go.Figure(go.Bar(x=bins,y=hist,))

        fig.update_layout(xaxis_title="Read length",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.html_L_hist, full_html=False, include_plotlyjs="cdn")

        # histogram of lengths with filtered out very long reads to facilitate histogram building at precision of 1bp
        series = in_data["lengths"].loc[lambda x: x <= 2500]

        nbins = max(series) - min(series)
        hist, bins = histogram(series,bins=nbins)
        fig = go.Figure(go.Bar(x=bins,y=hist,))

        fig.update_layout(xaxis_title="Read length",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)

        fig.write_html(output.html_L_hist_restrict, full_html=False, include_plotlyjs="cdn")

rule trimming:
    input:
        primers_path=primers_path,
        fastq = rules.filtQ.output
    output:
        fastq1="trimming/trim1.fastq",
        fastq2="trimming/trim2.fastq"
    shell:
        "ml load dorado/0.8.0 ; "
        "dorado trim --emit-fastq --primer-sequences {input.primers_path} {input.fastq} > {output.fastq1} ; "
        "dorado trim --emit-fastq --primer-sequences {primers_path} {output.fastq1} > {output.fastq2}"

rule filtL:
    input:
        rules.trimming.output.fastq2
    output:
        fastq="filtL/out.fastq"
    shell:
        "ml load nanofilt ; "
        "NanoFilt --length {min_len_ITS2} --maxlength {max_len_ITS2} {input} > {output}"

rule post_filtL_QC_nanoplot:
    input:
        rules.filtL.output
    output:
        pickle="filtL/nanoplot/NanoPlot-data.pickle"
    threads: 10
    shell:
        "ml load nanoplot ; "
        "NanoPlot -o filtL/nanoplot/ --fastq {input} --tsv_stats -t 10 --raw --store ; "

# rule post_filtL_QC_seqkit:
#     input:
#         rules.filtL.output
#     output:
#         txt="filtQ/seqkit_stats.txt"
#     threads: 1
#     shell:
#         "ml load seqkit ; "
#         "seqkit stats {input} > {output.txt} "

rule post_filtL_QC:
    input:
        pickle=rules.post_filtL_QC_nanoplot.output.pickle
    output:
        hist_L_html="filtL/hist_L.html",
        read_counts_table = "filtL/read_counts_table.pickle"
    run:
        from pickle import load
        from numpy import histogram
        import plotly.graph_objects as go
        import pandas
        from pickle import dump

        with open(input.pickle,"rb") as in_stream:
            in_data = load(in_stream)

        # calculate reads counts and export to pickle
        headers = ["", "All Q", f"Q > {minQ}"]
        out_list = []
        out_list.append((
            "reads in Total", len(in_data['lengths']), len(in_data['lengths'][(in_data['quals'] > minQ)])))
        out_list.append(("reads in range ITS2", len(
            in_data['lengths'][
                (in_data['lengths'] >= min_len_ITS2 ) & (in_data['lengths'] <= max_len_ITS2 )]), len(
            in_data['lengths'][
                (in_data['lengths'] >= min_len_ITS2 ) & (in_data['lengths'] <= max_len_ITS2 ) & (
                        in_data['quals'] > minQ)])))

        with open(output.read_counts_table,"wb") as out_stream:
            dump(out_list,out_stream)

        series = in_data["lengths"]

        nbins = max(series) - min(series)
        hist, bins = histogram(series,bins=nbins)
        fig = go.Figure(go.Bar(x=bins,y=hist,))

        fig.update_layout(xaxis_title="Read length",yaxis_title="Count",bargap=0)
        fig.update_traces(marker_line_width=0)
        fig.write_html(output.hist_L_html)

rule EMU:
    input:
        rules.filtL.output.fastq
    output:
        # log="EMU/out.log",
        abundance="EMU/out_rel-abundance.tsv",
        read_assignment_probs="EMU/out_read-assignment-distributions.tsv",
    conda:
        "custom-emu"
    threads: 20
    shell:
        # "conda info --envs > EMU/out.log 2>&1 ; "
        # "conda init ; conda activate /scratch/thdelcourt/metabarcoding/20250704/scripts/miniconda3/envs/custom-emu ; " #source /scratch/thdelcourt/metabarcoding/20250704/scripts/virtenv/EMU_ppl_py3.10/bin/activate ; "
        # filtering on all species with abundance > 0.0 because the default is 0.0001, which gets applied to the reads assignment file, leading to missed species. 
        "emu abundance --type map-ont --db {EMU_DB_path} "
        "--output-dir EMU  "
        "--keep-files --keep-counts --keep-read-assignments --output-unclassified --threads 20 --min-abundance 0.0 "
        "{input} > EMU/out.log 2>&1 ; "
        # "conda deactivate "

rule get_refs_for_consensus:
    input:
        rel_abundance_tsv="EMU/out_rel-abundance.tsv",
        sequences_dict = sequences_dict
    output:
        "refseqs.fasta"

    run:
        # from helper_scripts import check_to_barcode
        import pickle
        import csv

        # read in abundance data and check for each line if it needs to be realigned for consensus
        tax_ids = []
        with open(input.rel_abundance_tsv,"r") as in_stream:
            reader = csv.DictReader(in_stream,delimiter="\t")
            for record in reader:
                if record["tax_id"] not in ("unmapped", "mapped_unclassified"):  # skipping last 2 lines
                    # if check_to_barcode(record):
                    tax_ids.append(record["tax_id"])
        # print(tax_ids)

        with open(input.sequences_dict,"rb") as in_stream:
            in_dict_species_ids = pickle.load(in_stream)
        ref_seqs = []
        for tax_id in tax_ids:
            ref_seqs.extend(in_dict_species_ids[tax_id])
        # print(ref_seqs)

        # write a file with reference sequences taken from the DB of ref seqs
        out_data = ""
        for ref_seq in ref_seqs:
            with open(f"{ref_seqs_DB_path}/{ref_seq}.fasta","r") as in_stream:
                for line in in_stream.readlines():
                    out_data += line
            out_data += "\n"

        with open(f"{output}","w") as out_stream:
            out_stream.write(out_data)

rule generate_consensus:
    """
    Align fastq to references, then samtools formatting, then samtools consensus. Using minimap2 settings like EMU, and not removing secondary alignments to keep all info."""
    input:
        fasta=rules.get_refs_for_consensus.output,
        fastq=rules.filtL.output.fastq
    output:
        bam="consensus/alg_sort.bam",
        consensus_fasta="consensus/out_consensus_IUPAC.fa"
        # aligned_fasta="consensus/alg.fa"
    threads: 20
    shell:
        # exploration of data didn't reveal any QCFAIL or DUP reads. There were a few UNMAP reads, but they wouldn't end up in the pileup anyway. So the samtools --ff is essentially useless here.
        "ml load minimap2 ; "
        "ml load samtools ; "
        "minimap2 -ax 'map-ont' -t 10 -N  50 -p .9 -K 500M --secondary-seq {input.fasta} {input.fastq} | samtools view -Sb -@5 | samtools sort -@5 > {output.bam} ; "
        "samtools consensus --qual-calibration :r10.4_sup --homopoly-fix --homopoly-score 0.3 --low-MQ 5 --scale-MQ 1.5 "
        "--ff 'UNMAP,QCFAIL,DUP' -f fasta {output.bam} -A -o {output.consensus_fasta} ; "
        
        

rule count_aligned_reads_consensus:
    """
    Count number of total and aligned reads to all reference sequences.
    """
    input:
        bam=rules.generate_consensus.output.bam,
    output:
        stats_txt="consensus/samtools_stats_SN.txt",
        ref_counts="consensus/ref_counts.txt"
    threads: 1
    shell:
        "ml load samtools ; samtools index {input.bam} ; samtools stats -@ 5 {input.bam} | grep ^SN | cut -f 2- > {output.stats_txt} ; "  # calculate statistics and keep only SN (summary) metrics

        # iterates over sequences in the header and filters out primary alignments for each one of the references; reports the count of reads.
        # NOTE: filtering out the supplementary and secondary alignments is probably a bug as their annotation is (probably mostly) arbitrary. TO FIX: ensure to select all reads that provide correct information for the consensus generation step.
        # additionally, 0x40 referes to first in pair, which cannot be set in ONT, so it's useless here.
        # so don't filter on -F 0x40 -F 0x904
        'for seq in $(samtools view {input.bam} -H | grep "@SQ" | cut -f2 | cut -d":" -f2,3,4 --output-delimiter=":") ; do COUNT=$(samtools view -@ 5 {input.bam} $seq -c ) ; echo $seq $COUNT >> {output.ref_counts} ; done'
        #"samtools view -F 0x40 -F 0x904 -h -b consensus/alg_sort.bam > consensus/primary_alignments.bam ; " # keep only primary alignments
        #"samtools view primary.sam -H | grep '@SQ' | cut -f2 | cut -d':' -f2,3,4 - -output - delimiter=':' ; " # get all refseqs from header
        #"while read seq ; do echo $seq ; samtools view primary.bam $seq -c -F 0x40 -F 0x904 > consensus/refs/$seq.txt; done < consensus/ref_seqs.txt" # keep primary alignments and filter on each ref seq


rule collect_read_counts:
    input:
        input_reads_pickle=rules.plot_pre_filtQ_QC.output.read_counts_table,
        filtQ_reads_pickle=rules.plot_post_filtQ_QC.output.read_counts_table,
        filtL_reads_pickle=rules.post_filtL_QC.output.read_counts_table,
        EMU_rel_abundance_table=rules.EMU.output.abundance,
        consensus_aligned_reads_txt=rules.count_aligned_reads_consensus.output.stats_txt,
        consensus_aligned_reads_per_ref_txt=rules.count_aligned_reads_consensus.output.ref_counts
    output:
        txt="output/reads_counts_tables.txt",
    run:
        from pickle import load
        import tabulate
        import csv

        out_txt = ""
        with open(input.input_reads_pickle,"rb") as in_stream:
            in_data = load(in_stream)
            n_total_reads = in_data[0][1]
            n_total_ITS2_reads = in_data[1][1]
            # out_txt+=tabulate.tabulate(in_data,headers=headers)

        with open(input.filtQ_reads_pickle,"rb") as in_stream:
            in_data = load(in_stream)
            n_filtQ_reads = in_data[0][1]
            # out_txt += tabulate.tabulate(in_data,headers=headers)

        with open(input.filtL_reads_pickle,"rb") as in_stream:
            in_data = load(in_stream)
            n_filtL_reads = in_data[0][1]
            # out_txt += tabulate.tabulate(in_data,headers=headers)

        with open(input.EMU_rel_abundance_table,"r") as in_stream:
            reader = csv.DictReader(in_stream,delimiter="\t")
            EMU_est_counts=0

            for record in reader:
                if record["tax_id"] == "unmapped":
                    EMU_unmapped = int(float(record["estimated counts"]))
                elif record["tax_id"] == "mapped_unclassified":
                    EMU_mapped_unclassified = int(float(record["estimated counts"]))
                else:
                    EMU_est_counts+=float(record["estimated counts"])

        with open(input.consensus_aligned_reads_txt, 'r') as in_stream:
            in_data = in_stream.readlines()
            for line in in_data:
                line = line.strip().split("\t")
                if line[0]=="reads mapped:":
                    consensus_mapped=int(line[1])
                if line[0] == "reads unmapped:":
                    consensus_unmapped = int(line[1])

        out_list = []
        headers=["Data", "Reads counts"]
        out_list.append(["Input", n_total_reads])
        out_list.append(["Input ITS2 size", n_total_ITS2_reads])
        out_list.append(["FiltQ", n_filtQ_reads])
        out_list.append(["FiltL", n_filtL_reads])
        out_list.append(["EMU estimate", int(EMU_est_counts)])
        out_list.append(["EMU unmapped", EMU_unmapped])
        out_list.append(["EMU mapped unclassified", EMU_mapped_unclassified])
        out_list.append(["Consensus mapped", consensus_mapped])
        out_list.append(["Consensus unmapped", consensus_unmapped])

        with open(output.txt, 'w') as out_stream:
            out_stream.write(tabulate.tabulate(out_list, headers=headers))


rule viralMSA:
    """Resorting to the tmp_dir and mv hack because ViralMSA errors if the output dir exists,
    and snakemake creates it automatically."""
    input:
        consensus_fa=rules.generate_consensus.output.consensus_fasta,
        path_to_Hoodia_SNP_barcode_ref=path_to_Hoodia_SNP_barcode_ref
    output:
        dir=directory("ViralMSA"),
        aln_fasta="ViralMSA/out_consensus_IUPAC.fa.aln"
    shell:
        "ml load minimap2 ; "
        "{path_to_ViralMSA} -s {input.consensus_fa} -r {input.path_to_Hoodia_SNP_barcode_ref} -o {output.dir}/tmp ; "
        "mv {output.dir}/tmp/* {output.dir}"

rule extract_barcode:
    input:
        rel_abundance_tsv = "EMU/out_rel-abundance.tsv",
        fasta=rules.viralMSA.output.aln_fasta
    output:
        barcode_results="SNP_barcode/output/SNP_barcode_results.pckl"
    run:
        from Bio import SeqIO
        from pickle import dump
        from yaml import safe_load
        with open(path_to_ITS2_barcode_positions, "r") as in_stream:
            ITS2_barcode_positions=safe_load(in_stream)

        from helper_scripts import check_to_barcode

        import csv

        # read in abundance data and check for each line if it needs to be realigned for consensus
        tax_ids = set()
        with open(input.rel_abundance_tsv,"r") as in_stream:
            reader = csv.DictReader(in_stream,delimiter="\t")
            for record in reader:
                if record["tax_id"] not in ("unmapped", "mapped_unclassified"):  # skipping last 2 lines
                    if check_to_barcode(record):
                        tax_ids.add(record["tax_id"])


        aligned_seqs = SeqIO.to_dict(SeqIO.parse(input.fasta,"fasta"))
        SNP_barcode_results = []
        # results_hoodias = set()
        # results_non_hoodias = set()
        for alg_key in aligned_seqs:
            if alg_key.split(":")[0] in tax_ids:
                record = aligned_seqs[alg_key]
                seq = record.seq
                # print(seq)
                indexes = list(ITS2_barcode_positions.keys())
                barcode = tuple(seq[pos-1].upper() == ITS2_barcode_positions[pos] for pos in ITS2_barcode_positions.keys())
                seq_name = record.description
                SNP_barcode_results.append((seq_name, barcode))
                # if all(barcode):
                #     results_hoodias.add(seq_name)
                # else:
                #     results_non_hoodias.add(seq_name)
            # print(SNP_barcode_results)

            with open(output.barcode_results, "wb") as out_stream:
                dump(SNP_barcode_results, out_stream)

rule get_consensus_alignments_identity:
    """
    Align and analyse the consensus generated to calculate identity to the expected ref sequence for specific species.  
    """
    input:
        ref_fasta = rules.get_refs_for_consensus.output,
        consensus_fasta = rules.generate_consensus.output.consensus_fasta,
    output:
        txt = "consensus_identity/out.tsv",
        pickle = "consensus_identity/out.pickle"
    run:
        from Bio import Align # substitution_matrices
        from Bio import SeqIO
        from pickle import dump
        from helper_scripts import align_2_seqs, count_identical_bases, calc_id
        # for debugging
        import traceback


        # define aligner
        # strives to be balanced.
        # set as needle emboss: matrix EDNAFULL (=NUC.4.4)
        # + scores for gaps
        # But reduce end gaps scores compared to gaps to help end-gaps to occur.
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = Align.substitution_matrices.load("NUC.4.4")
        # #use global alignment for ease of parsing.
        aligner.mode = "global"
        # # scores aim at avoiding internal gaps, but allowing end gaps due to different amplicon construction.
        aligner.internal_open_gap_score = -10
        aligner.internal_extend_gap_score = -0.5
        aligner.end_open_gap_score = -5
        aligner.end_extend_gap_score = -0.1

        # debugging error-less bug / hanging process
        try:
            print("Aligner definition:")
            print(aligner)
            # for i in aligner.substitution_matrix:
            #     print(i)
            consensus_sequences = SeqIO.to_dict(SeqIO.parse(input.consensus_fasta,"fasta"))
            ref_seqs = SeqIO.to_dict(SeqIO.parse(str(input.ref_fasta),"fasta"))
    
            # Filter input sequences to keep only consensus sequences that are of good enough quality (eg not full "NNNN") - probably not that usefu now that minimap outputs all alignment nucleotides
            cons_seq_to_use = [cons_seq for cons_seq in consensus_sequences if
                               sum([nt not in ("A", "T", "C", "G") for nt in consensus_sequences[cons_seq].seq]) != len(
                                   consensus_sequences[cons_seq].seq)]
    
            # create alignment pairs list
            input_list = []
            for cons_seq in cons_seq_to_use:
                target = ref_seqs[cons_seq]
                query = consensus_sequences[cons_seq]
                identifier = target.id  # no need for the query.id as it's the same sequence name
                                        # as we're aligning the consensus of 1 ref to the ref
                input_list.append((aligner, identifier, target.seq.upper(), query.seq.upper()))
            
            print(f"Input list length: {len(input_list)}")
            # align sequence pairs; not in parallel as it is not expected that there will be that many sequences to align. Max 1-2 seconds.

            results_list = [align_2_seqs(input_data) for input_data in input_list]

            print(f"Results list length: {len(input_list)}")
            
            # group results per species in output dictionary
            dict_species = dict()
            for result in results_list:
                name = result[0]
                species_id = name.split(":")[0]
                if species_id not in dict_species:
                    dict_species[species_id] = []
                dict_species[species_id].append(result)
            print(f"dict_species length: {len(dict_species)}")
            
            # calculate alignment scores of all alignments per species for further filtering and reporting
            for species in dict_species:
                # best_score = 10000000   # arbitrarily big as there won't be alignments with that many differences in a =250bp sequence.
                # best_alignment = None
                tmp_list = []
                for identifier, alignment in dict_species[species]:
                    identity = count_identical_bases(alignment[0], alignment[1])
                    clustalID = calc_id(alignment[0], alignment[1])
                    diffs = len(alignment[0]) - identity    # take any sequence as they are of the same length after aligning
                    tmp_list.append((identifier, alignment, identity, clustalID, diffs))
                dict_species[species] = sorted(tmp_list, key=lambda record: record[4], reverse=False) # add to dict as sorted list by diffs value (smaller = better)
            
            # export pickle dictionary
            with open(output.pickle, "wb") as out_stream:
                dump(dict_species, out_stream)
    
            # export full alignment text file
            with open(output.txt, "w") as out_stream:
                for species in dict_species:
                    out_stream.write(f"\n#######\n{species}\n#######\n")
                    for identifier, alignment, identity, clustalID, diffs in dict_species[species]:
                        out_stream.write(identifier)
                        out_stream.write(f"identity: {identity}    differences: {diffs}    clustalID: {clustalID}\n")
                        out_stream.write(alignment.format(fmt=""))
                        out_stream.write("\n")
        except Exception:
            print(traceback.format_exc())
            
rule group_identical_species_and_combine_outputs:
    """
    Group species with identical sequence and add their abundances together.
    Combine abundance and barcode output.
    Generates a table of combined outputs and grouped records based on the EMU relative abundance table.
    """
    input:
        abundance=rules.EMU.output.abundance,
        taxonomy_pickle=taxonomy_pickle,
        SNP_barcode_file_path=rules.extract_barcode.output.barcode_results,
        realignment_consensus_file_path=rules.get_consensus_alignments_identity.output.txt,
        consensus_alignment_pickle = rules.get_consensus_alignments_identity.output.pickle,

    output:
        tsv="output/results.tsv",
    run:
        from pickle import load
        # from helper_scripts import get_group_from_tax_id
        # from helper_scripts import get_tax_ids_from_group
        from helper_scripts import group_identical_species_and_combine_outputs
        from  yaml import safe_load

        with open(config["path_to_identical_species_DB"],'r') as in_stream:
            dict_groups_values = safe_load(in_stream)

        # create dict_species_to_group by 'inverting' the dict_groups_values
        dict_species_to_group = dict()
        for group in dict_groups_values:
            for species in dict_groups_values[group]:
                dict_species_to_group[species] = group

        with open(input.SNP_barcode_file_path, "rb") as in_stream:
            SNP_barcode_data = load(in_stream)

        # with open(input.consensus_alignment_pickle, "rb")
        # print(output.tsv)
        output_data = group_identical_species_and_combine_outputs(abundance_file_path=input.abundance, SNP_barcode_data=SNP_barcode_data,
            taxonomy_pickle_path=taxonomy_pickle, output_path=output.tsv,
            dict_species_to_group=dict_species_to_group, Hoodia_ref=Hoodia_ref,
            dict_groups_values=dict_groups_values, consensus_alignment_pickle_path = input.consensus_alignment_pickle)

rule tabular_ref_counts_and_cons_diffs:
    """
    Generate tsv table of counts of aligned reads to reference sequences and differences from consensus for each ref sequence.
    columns: species id  |  species name  |  reference seq id  |  aligned read counts  |  diff between consensus and reference
    rows: 1 per ref sequences of species in EMU output (may be multiple seq per species)
    """
    input:
        consensus_alignment_pickle = rules.get_consensus_alignments_identity.output.pickle,
        taxonomy_pickle = taxonomy_pickle,
        ref_counts_txt = rules.count_aligned_reads_consensus.output.ref_counts,
    output:
        table="output/tabular_ref_counts_and_diffs.tsv",
        
    run:
        from csv import writer
        from pickle import load
        
        # load taxonomy for species names
        with open(input.taxonomy_pickle,"rb") as in_stream:
            dict_taxonomy = load(in_stream)

        # load consensus alignment pickle for consensus alignment infos 
        with open(input.consensus_alignment_pickle,'rb') as in_stream:
            dict_consensus_alignment = load(in_stream)

        # load ref counts txt for ref counts info and make it a dict of counts per ref seq name
        with open(input.ref_counts_txt,'r') as in_stream:
            dict_ref_counts = {record[0]: int(record[1]) for record in
                               [line.strip("\n").split(" ") for line in in_stream.readlines()]}

        output_data = []
        for tax_id in dict_consensus_alignment:

            species_name = dict_taxonomy[tax_id]

            for consensus_alignment_record in dict_consensus_alignment[tax_id]:
                ref_seq_name = consensus_alignment_record[0]
                ref_cons_al_len = consensus_alignment_record[2]
                ref_cons_al_id = round(consensus_alignment_record[3] * 100,2)
                ref_cons_al_diffs = consensus_alignment_record[4]
                ref_seq_counts = dict_ref_counts[ref_seq_name]

                output_data.append([tax_id, species_name, ref_seq_name, ref_seq_counts, ref_cons_al_diffs,
                                    ref_cons_al_len, ref_cons_al_id])

        # write tsv output file
        header = ["tax_id", "species_name", "ref_seq_name", "ref_seq_counts", "ref_cons_al_diffs", "ref_cons_al_len",
                  "ref_cons_al_id"]
        with open(output.table,'w') as csvfile:
            csvwriter = writer(csvfile,delimiter='\t')
            csvwriter.writerow(header)
            csvwriter.writerows(output_data)
        
        
rule read_assignment_probabilities:
    """
    Get mean and median values of read assignment probabilities per species. 
    """
    input:
        read_assignment_probs=rules.EMU.output.read_assignment_probs,
        taxonomy_pickle= taxonomy_pickle,
    output:
        read_assignment_stats="output/read_assignment_stats.tsv",
    run:
        from pandas import DataFrame
        from pandas import read_csv
        from pickle import load

        with open(taxonomy_pickle,"rb") as in_stream:
            dict_taxonomy = load(in_stream)
        
        df_reads = read_csv(input.read_assignment_probs,sep="\t", index_col=0)
        
        dict_sp_names = dict()
        for tax_id in df_reads:
            try:
                dict_sp_names[tax_id] = dict_taxonomy[tax_id]
            except:
                dict_sp_names[tax_id] = "NA"
                
        df_stats = DataFrame.from_dict(dict_sp_names,orient='index',columns=["species_name"])
        df_stats["median"] = df_reads.median()
        df_stats["mean"] = df_reads.mean()
        df_stats = df_stats.round(decimals=7)
        df_stats.to_csv(path_or_buf=output.read_assignment_stats,sep='\t')
        
    
rule tabular_QCs:
    """
    Collate various metrics into tabular output format: read counts before and after filtering and alignment."""
    input:
        prefiltQ_pickle=rules.pre_filtQ_QC_nanoplot.output.pickle,
        post_filtQ_pickle=rules.post_filtL_QC_nanoplot.output.pickle,

    output:
        table="output/tabular_QC.tsv",
    run:
        from pickle import load
        from tabulate import tabulate

        with open(input.prefiltQ_pickle,"rb") as in_stream:
            in_data = load(in_stream)
        headers = ["", "All Q", f"Q > {minQ}"]
        out_list = []
        out_list.append(("reads in Total", len(in_data['lengths']),
                         len(in_data['lengths'][(in_data['quals'] > minQ)])))
        out_list.append(("reads in range ITS2", len(
            in_data['lengths'][(in_data['lengths'] >= min_len_ITS2) & (in_data['lengths'] <= max_len_ITS2)]), len(
            in_data['lengths'][(in_data['lengths'] >= min_len_ITS2) & (in_data['lengths'] <= max_len_ITS2) & (
                        in_data['quals'] > minQ)])))
        with open(output.table, "w") as out_stream:
            out_stream.write("prefiltQ\n")
            out_stream.write(tabulate(out_list,headers=headers))
            out_stream.write("\n")


        with open(input.post_filtQ_pickle,"rb") as in_stream:
            in_data = load(in_stream)
            headers = ["", "All Q", f"Q > {minQ}"]
            out_list = []
            out_list.append(("reads in Total", len(in_data['lengths']),
                             len(in_data['lengths'][(in_data['quals'] > minQ)])))
            out_list.append(("reads in range ITS2", len(
                in_data['lengths'][(in_data['lengths'] >= min_len_ITS2) & (in_data['lengths'] <= max_len_ITS2)]), len(
                in_data['lengths'][(in_data['lengths'] >= min_len_ITS2) & (in_data['lengths'] <= max_len_ITS2) & (
                            in_data['quals'] > minQ)])))

        with open(output.table,"a") as out_stream:
            out_stream.write("\npost_filtQ\n")
            out_stream.write(tabulate(out_list,headers=headers))
            out_stream.write("\n")

rule filter_EMU_output:
    """
    Filters species identified by EMU based on 4 filters.
    """
    
    input:
        emu_results_path=rules.EMU.output.abundance,
        reads_consensus_file_path=rules.tabular_ref_counts_and_cons_diffs.output.table,
        reads_assignment_file_path=rules.read_assignment_probabilities.output.read_assignment_stats,
    output:
        filtered_table="output/filtered_species.tsv",
        full_table="output/full_species.tsv",
    run:
        import pandas 
        
        def load_and_combine_data(EMU_data_path=None, reads_consensus_file_path = None, reads_assignment_file_path=None):
            """
            Load data from out_rel-abundance.tsv, tabular_ref_counts_and_diffs.tsv and read_assignment_stats.tsv.
            Extract the data, clean it if necessary (remove unclassified and unmapped reads), select the consensus data for the ref with max alignments counts 
            (or lowest differences if ties, or first in line if both metrics have ties) and combine the data per EMU output record (= identified species).
            
            :param EMU_data_path: 
            :param reads_consensus_file_path: 
            :param reads_assignment_file_path: 
            :return: merged dataframe
            """
            # load EMU output data:
            EMU_data_path = EMU_data_path 
            # f"{directory}/EMU/out_rel-abundance.tsv"
            df_EMU_data = pandas.read_csv(EMU_data_path,sep="\t")
            df_EMU_data_map_class = df_EMU_data[~df_EMU_data["tax_id"].isin(["unmapped", "mapped_unclassified"])]
            
            df_EMU_data_map_class = df_EMU_data_map_class.astype({"tax_id": "string"})

            # Load reads assignment metrics
            reads_consensus_file_path = reads_consensus_file_path
            # f"{directory}/output/tabular_ref_counts_and_diffs.tsv"
            df_reads_cons = pandas.read_csv(reads_consensus_file_path,sep="\t")
            df_reads_cons = df_reads_cons.astype({"tax_id": "string"})

            # select the ref seq with highest reads counts only. ==> tax_id is now unique in df
            idx = df_reads_cons.groupby(['tax_id'])['ref_seq_counts'].transform("max") == df_reads_cons[
                'ref_seq_counts']
            df_reads_cons_max_ref_counts = df_reads_cons[idx]

            # in cases where more than one ref remains (low read counts, it can happen), also filter on differences,and then take first one of group. happened in 2 cases
            idx = df_reads_cons_max_ref_counts.groupby(['tax_id'])['ref_cons_al_diffs'].transform("min") == \
                  df_reads_cons_max_ref_counts['ref_cons_al_diffs']
            df_reads_cons_max_ref_counts = df_reads_cons_max_ref_counts[idx]
            df_reads_cons_max_ref_counts = df_reads_cons_max_ref_counts.groupby(['tax_id']).first()

            #Load reads assignment data
            reads_assignment_file_path = reads_assignment_file_path
            # f"{directory}/output/read_assignment_stats.tsv"
            df_reads_assgnt = pandas.read_csv(reads_assignment_file_path,sep="\t",header=0,names=["tax_id",
                                                                                                  "species_name",
                                                                                                  "median", "mean"])
            df_reads_assgnt = df_reads_assgnt.astype({"tax_id": "string"})

            # merge 3 dataframes using the tax_id as key, as it should be unique within each sample for these 3 dataframes.
            combined_df = pandas.merge(df_EMU_data_map_class,df_reads_cons_max_ref_counts.drop(columns=[
                "species_name"]),how="outer",on="tax_id")

            combined_df = pandas.merge(combined_df,df_reads_assgnt.drop(columns=[
                "species_name"]),how="outer",on="tax_id")

            return combined_df
        
        def add_filter_column(data=None, abundance_thr=None, read_assign_med_thr=None, ref_cov_thr=None, consensus_diffs_thr=None):
            """
            Add a column with filtered out information based on 4 threshold values.
            SHOULD BE REWORKED TO BE MORE PANDA-ISH INSTEAD OF ITERATING PER ROW
            :param data: dataframe
            :param abundance_thr: abundance threshold
            :param read_assign_med_thr: read assignment threshold
            :param ref_cov_thr: reference coverage threshold
            :param consensus_diffs_thr: consensus differences threshold
            :return: dataframe with added column to filter on
            """
            
            data["filter_out"] = False
            for index in range(len(data)):

                a_filter_fails = False
    
                # filter on abundance
                abundance = float(data.loc[index, "abundance"])
                if pandas.isna(abundance):
                    abundance = 0
                else:
                    abundance = float(abundance)
    
                if abundance < abundance_thr:
                    a_filter_fails = True
    
                # filter on read assignment median
                read_assignment_median = data.loc[index, "median"]
                if pandas.isna(read_assignment_median):
                    read_assignment_median = 0
                else:
                    read_assignment_median = float(read_assignment_median)
                
                if read_assignment_median < read_assign_med_thr:
                    a_filter_fails = True
    
                # filter on ref coverage
                ref_cov = data.loc[index, "ref_seq_counts"]
                if pandas.isna(ref_cov):
                    ref_cov = 0
                else:
                    ref_cov = int(ref_cov)
                
                if ref_cov < ref_cov_thr:
                    a_filter_fails = True
    
                # filter on consensus differences
                cons_diff = data.loc[index, "ref_cons_al_diffs"]
                if pandas.isna(cons_diff):
                    cons_diff = 500  # expected max based on 250 sequences alignments; well beyond anything useful anyway
                else:
                    cons_diff = int(cons_diff)
    
                if cons_diff > consensus_diffs_thr:
                    a_filter_fails = True
    
                if a_filter_fails:
                    data.loc[index, "filter_out"] = True
            
            return data
        
        combined_df = load_and_combine_data(EMU_data_path=input.emu_results_path, reads_consensus_file_path = input.reads_consensus_file_path, reads_assignment_file_path=input.reads_assignment_file_path)
        df_to_filter = add_filter_column(data=combined_df, abundance_thr=min_abundance, read_assign_med_thr=min_median_rap, ref_cov_thr=min_real_depth, consensus_diffs_thr=max_diff_cons_ref)
        
        filtered_df = df_to_filter[df_to_filter["filter_out"] == False]

        filtered_df.to_csv(path_or_buf=output.filtered_table,sep='\t',index=True,decimal=".")
        df_to_filter.to_csv(path_or_buf=output.full_table,sep='\t',index=True,decimal=".")
        
        
rule html_report:
    """Writes a basic report with output results. This creates a self-contained html (quite heavy).
    Ultimately, maybe better to use the snakemake solution using jinja2 file definitions.
    """
    input:
        table_QC = rules.tabular_QCs.output.table,
        # table_abundance = rules.group_identical_species_and_combine_outputs.output.tsv,
        table_abundance= rules.filter_EMU_output.output.filtered_table,

        pre_html_Q_hist=rules.plot_pre_filtQ_QC.output.html_Q_hist,
        pre_html_Q_cum_hist=rules.plot_pre_filtQ_QC.output.html_Q_cum_hist,
        pre_html_L_hist=rules.plot_pre_filtQ_QC.output.html_L_hist,
        pre_html_L_hist_restrict=rules.plot_pre_filtQ_QC.output.html_L_hist_restrict,
        # pre_txt_seqkit=rules.pre_filtQ_QC_seqkit.output.txt,

        post_filtQ_html_Q_hist=rules.plot_post_filtQ_QC.output.html_Q_hist,
        post_filtQ_html_Q_cum_hist=rules.plot_post_filtQ_QC.output.html_Q_cum_hist,
        post_filtQ_html_L_hist=rules.plot_post_filtQ_QC.output.html_L_hist,
        post_filtQ_html_L_hist_restrict=rules.plot_post_filtQ_QC.output.html_L_hist_restrict,
        # post_filtQ_txt_seqkit=rules.post_filtQ_QC_seqkit.output.txt,

        post_filtL_hist_L_html=rules.post_filtL_QC.output.hist_L_html,
        # psot_filtL_seqkit=rules.post_filtL_QC_seqkit.output.txt,

        reads_counts_table=rules.collect_read_counts.output.txt,

    output:
        html="output/report.html"
    run:
        from time import strftime

        with open(input.table_QC,"r") as in_stream:
            table_QC = in_stream.read()

        with open(input.reads_counts_table, "r") as in_stream:
            table_reads_counts = in_stream.read()

        # with open(input.pre_txt_seqkit, "r") as in_Stream:
        #     table_pre_txt_seqkit = in_stream.read()

        with open(input.table_abundance,"r") as in_stream:
            table_abundance = in_stream.read()

        with open(input.pre_html_Q_hist,"r") as in_stream:
            pre_html_Q_hist = in_stream.read()

        with open(input.pre_html_Q_cum_hist,"r") as in_stream:
            pre_html_Q_cum_hist = in_stream.read()

        with open(input.pre_html_L_hist,"r") as in_stream:
            pre_html_L_hist = in_stream.read()

        with open(input.pre_html_L_hist_restrict,"r") as in_stream:
            pre_html_L_hist_restrict = in_stream.read()

        with open(input.post_filtQ_html_Q_hist,"r") as in_stream:
            post_filtQ_html_Q_hist = in_stream.read()

        with open(input.post_filtQ_html_Q_cum_hist,"r") as in_stream:
            post_filtQ_html_Q_cum_hist = in_stream.read()

        with open(input.post_filtQ_html_L_hist,"r") as in_stream:
            post_filtQ_html_L_hist = in_stream.read()

        with open(input.post_filtQ_html_L_hist_restrict,"r") as in_stream:
            post_filtQ_html_L_hist_restrict = in_stream.read()

        with open(input.post_filtL_hist_L_html,"r") as in_stream:
            post_filtL_hist_L_html = in_stream.read()

        generation_time = strftime("%Y-%m-%d %H:%M:%S")

        html_string = ('''
        <html>
            <head>
                <title>Analysis report for '''+ sample_name + '''</title>
                <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                <style>
                    body { margin:0 100; background:whitesmoke; }
                </style>
            </head>
            <body>
                <h1>Analysis report for '''+ sample_name + '''</h1>
                <p>Generated on ''' + generation_time + ''' </p>
                <h2>Results:</h2>
                <h3>Abundance and barcode info</h3>
                <pre>''' + table_abundance + '''</pre>
                <h2>QC:</h2>
                <h3>Read counts</h3>
                <pre>''' + table_reads_counts + '''</pre>
                <h3>Read counts pre and post trimming</h3>
                <pre>'''+ table_QC +'''</pre>
                <h3>Pre filtering Q</h3>
                <div style="height: 600px; width:1200px;" >
                ''' + pre_html_Q_hist + '''
                </div>
                <div style="height: 600px; width:1200px;" >
                ''' + pre_html_Q_cum_hist + '''
                </div>
                <div style="height: 600px; width:1200px;" >
                ''' + pre_html_L_hist + '''
                </div>
                <div style="height: 600px; width:1200px;" >
                ''' + pre_html_L_hist_restrict + '''
                </div>
                
                <h3>Post filtering Q</h3>
                <p>Before trimming!</p>
                <div style="height: 600px; width:1200px;" >
                ''' + post_filtQ_html_Q_hist + '''
                </div>
                <div style="height: 600px; width:1200px;" >
                ''' + post_filtQ_html_Q_cum_hist + '''
                </div>
                <div style="height: 600px; width:1200px;" >
                ''' + post_filtQ_html_L_hist + '''
                </div>
                <div style="height: 600px; width:1200px;" >
                ''' + post_filtQ_html_L_hist_restrict + '''
                </div>
                
                <h3>Post filtering Length</h3>
                <div style="height: 600px; width:1200px;" >
                ''' + post_filtL_hist_L_html + '''
                </div>
            </body>
        </html>''')

        with open(output.html,'w') as out_stream:
            out_stream.write(html_string)




rule all:
    input:
        # QC
        rules.pre_filtQ_QC_nanoplot.output.pickle,
        rules.plot_pre_filtQ_QC.output.html_Q_hist,
        rules.plot_pre_filtQ_QC.output.html_Q_cum_hist,
        rules.plot_pre_filtQ_QC.output.html_L_hist,
        rules.plot_pre_filtQ_QC.output.html_L_hist_restrict,
        # rules.post_filtQ_QC_seqkit.output.txt,
        rules.plot_post_filtQ_QC.output.html_Q_hist,
        rules.plot_post_filtQ_QC.output.html_Q_cum_hist,
        rules.plot_post_filtQ_QC.output.html_L_hist,
        rules.plot_post_filtQ_QC.output.html_L_hist_restrict,
        rules.post_filtL_QC.output.hist_L_html,
        # rules.post_filtL_QC_seqkit.output.txt,
        rules.tabular_QCs.output.table,
        rules.collect_read_counts.output.txt,
        rules.tabular_ref_counts_and_cons_diffs.output.table,
        rules.read_assignment_probabilities.output.read_assignment_stats,
        rules.filter_EMU_output.output.filtered_table,
        rules.filter_EMU_output.output.full_table,
        
        # final results
        rules.group_identical_species_and_combine_outputs.output,
        rules.html_report.output
    default_target: True


