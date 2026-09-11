def check_to_barcode(record=None, decision_base="Apocynaceae"):
    """
    Check if the provided EMU output record is from a species that warrants SNP-barcoding.
    Currently, if the species is from Apocynaceae ==> yes otherwise, no.
    Could also be based on the best identity between any sequences of that species and any Hoodia or Apocynaceae (pre-set list of species).
    record: dict
    decision_base: string [Apocynaceae, species_list]
    Record of the format (tab-separated):
    {'tax_id': '41737', 'abundance': 0.055081812721262245, 'kingdom': 'Viridiplantae', 'phylum': 'Streptophyta',
    'class': 'Magnoliopsida', 'order': 'Gentianales', 'family': 'Apocynaceae', 'genus': 'Hoodia',
    'species': 'Hoodia_officinalis', 'estimated counts': 253.81699301957642}
    """
    if decision_base == "Apocynaceae":
        if record["family"] == "Apocynaceae":
            return True
        else:
            return False

# def get_species_to_group(in_tax_id, dict_species_to_group=dict_species_to_group, dict_groups_values=dict_groups_values):
#     if in_tax_id in dict_species_to_group:
#         return dict_groups_values[dict_species_to_group[in_tax_id]]
#     else:
#         return [in_tax_id,]

def get_group_from_tax_id(in_tax_id, dict_species_to_group=None):
    if in_tax_id in dict_species_to_group:
        return dict_species_to_group[in_tax_id]
    else:
        return "-1"

def get_tax_ids_from_group(in_group, dict_groups_values=None):
    return dict_groups_values[in_group]


def group_identical_species_and_combine_outputs(abundance_file_path=None, SNP_barcode_data=None,
                                                taxonomy_pickle_path=None, output_path=None,
                                                dict_species_to_group=None, Hoodia_ref=None,
                                                dict_groups_values=None, consensus_alignment_pickle_path=None):
    """
    Group species with identical sequence and add their abundances together.
    Combine abundance and barcode output.
    Generates a table of combined outputs and grouped records based on the EMU relative abundance table.
    outputs a tsv table
    """

    import csv
    from tabulate import tabulate
    from pickle import load

    ### Import results data ###

    # read in EMU abundance data
    rel_abundance_results = []
    with open(abundance_file_path, "r", ) as in_stream:
        reader = csv.DictReader(in_stream, delimiter="\t")
        for line in reader:
            if line["tax_id"] not in ("unmapped", "mapped_unclassified"):  # skip last 2 lines
                line["abundance"] = float(line["abundance"])
                line["estimated counts"] = float(line["estimated counts"])
                rel_abundance_results.append(line)
    with open(taxonomy_pickle_path, "rb") as in_stream:
        dict_taxonomy = load(in_stream)
    with open(consensus_alignment_pickle_path, 'rb') as in_stream:
        dict_consensus_alignment = load(in_stream)
    print(dict_consensus_alignment)

    # load tax_ids and abundance data into groups of 'identical' species into species_groups_data
    # non-grouped species are added to group "-1" and treated differently down the line
    species_groups_data = dict()
    for line in rel_abundance_results:
        tax_id = line["tax_id"]
        tax_group = get_group_from_tax_id(tax_id, dict_species_to_group=dict_species_to_group)
        if tax_group not in species_groups_data:
            species_groups_data[tax_group] = {"abdce": [], "barcode": [], "consensus_diffs": []}
        species_groups_data[tax_group]["abdce"].append(line)

    # load barcoding results and append them to a list in the correct taxonomical group in species_groups_data
    for line in SNP_barcode_data:
        tax_id = line[0].split(":")[0]
        if tax_id != Hoodia_ref:  # avoid the ref which is the "rucschii" name -
                                  # Probably not useful as the ref isn't in the barcode data anymore at this point.
                                  # This could cause a bug if the names aren't modified in the DB.
                                  # But given that EMU uses a modified DB, it's fine.
            tax_group = get_group_from_tax_id(tax_id, dict_species_to_group=dict_species_to_group)
            species_groups_data[tax_group]["barcode"].append(line)

    # load consensus-based identity info (n diffs => [5] in list) into species grouped dictionary
    # (only the best alignment - aka the [0] in list)
    for tax_id in dict_consensus_alignment:
        # keep only the best alignment info per species
        cons_diff_val = min([i[4] for i in dict_consensus_alignment[tax_id]])
        tax_group = get_group_from_tax_id(tax_id, dict_species_to_group=dict_species_to_group)
        species_groups_data[tax_group]["consensus_diffs"].append((tax_id, cons_diff_val))


    ### Format and combine inputs into single output ###
    output_data = []

    for group in species_groups_data:

        # this way of dealing with grouped and non-grouped species isn't clean. To refactor at some point.

        if group == "-1":   # for non-grouped species; create new dict to reconcile barcode infos and abundance infos with tax_id
                            # and print without group id (-1)

            ungrouped_species_dict = dict()
            for result in species_groups_data[group]['abdce']:
                ungrouped_species_dict[result["tax_id"]] = [result, [], 10000]

            for seq_name, barcode_data in species_groups_data[group]["barcode"]:
                tax_id = seq_name.split(":")[0]
                ungrouped_species_dict[tax_id][1].append((seq_name, barcode_data))

            for tax_id, consensus_diff in species_groups_data[group]["consensus_diffs"]:
                ungrouped_species_dict[tax_id][2] = consensus_diff

            for tax_id in ungrouped_species_dict:

                pct = ungrouped_species_dict[tax_id][0]["abundance"]
                species_name = dict_taxonomy[tax_id]
                cons_diff = ungrouped_species_dict[tax_id][2]

                if ungrouped_species_dict[tax_id][1]==[]: # if there is no barcode data for this species, generate output without barcode

                    output_data.append(["", round(pct, 5), species_name, "", "", cons_diff])

                else:  # if there is barcode data for this species, generate output with barcode infos (may be multiple lines)

                    first_in_tax_id = True  # to only print species info on first line (first barcode result) to not repeat info

                    for seq_name, barcode_data in ungrouped_species_dict[tax_id][1]:
                        if first_in_tax_id:
                            first_in_tax_id = False
                            output_data.append(["", round(pct, 5), species_name, seq_name,
                                                ','.join(['1' if i is True  else '0' for i in barcode_data ]),
                                                cons_diff])
                        else:
                            output_data.append(["", '"', '"',
                                                seq_name, ','.join(['1' if i is True  else '0' for i in barcode_data ]),
                                                '"'])

        else:  # for grouped species

            pct = sum([line["abundance"] for line in species_groups_data[group]["abdce"]])
            group_species_names = [dict_taxonomy[tax_id] for tax_id in
                                   get_tax_ids_from_group(group, dict_groups_values=dict_groups_values)]

            consensus_diff = min([consensus_diff for tax_id, consensus_diff in species_groups_data[group]["consensus_diffs"]])

            if species_groups_data[group]["barcode"] == []:  # if there is no barcode data for this species, generate output without barcode info

                output_data.append([group, round(pct, 5), ', '.join(group_species_names), "", "", consensus_diff])

            else:
                first_in_group = True  # to only print group and species info once per group to avoid repetitions

                for seq_name, barcode_data in species_groups_data[group]["barcode"]:

                    if first_in_group:  # add full species and group info if first line (barcode) of group
                        output_data.append([group, round(pct, 5), ', '.join(group_species_names), seq_name,
                                            ','.join(['1' if i is True  else '0' for i in barcode_data ]),
                                            consensus_diff])
                        first_in_group = False
                    else:
                        output_data.append(['"', '"', '"', seq_name,
                                            ','.join(['1' if i is True  else '0' for i in barcode_data ]),
                                            '"'])

    header = ["Group id", "Gp species abundance sum", "Group species names", "Sequence name", "Barcode", "Consensus_diff"]
    with open(output_path, "w") as out_stream:
        out_stream.write(tabulate(output_data, headers=header, floatfmt=".5f"))

def align_2_seqs(args):
    """ALign 2 sequences with smith-waterman algorithm implemented in BioPython.
    Input list formatted as:
    [[aligner, target_name, query_name, target_seq, query_seq], []...]"""
    aligner, identifier, target, query = args[:]
    try:
        # store 'best' alignment (first of iterator; if multiple => same score so it doesn't matter much which one we take) and the score of the alignment
        aligts = aligner.align(target,query)
    except Exception:
        print("error: ")
        print(f"{identifier}")
        # traceback.print_exc()
        raise Exception
    return identifier, aligts[0], #aligts.score


def calc_id(seq1, seq2):
    """
    Calculate identity as clustal W:
        n ids of those that are not indel (big assumption that will bias results)
    seq1 and seq2 must be of same length.
    """
    if len(seq1) == len(seq2):
        n_ids = sum([1 for i in range(len(seq1)) if seq1[i] == seq2[i]])
        n_not_gaps = sum([1 for i in range(len(seq1)) if seq1[i] != "-" and seq2[i] != "-"])
        if n_not_gaps > 0:
            prop_identity = n_ids / n_not_gaps
        else:
            prop_identity = 0
    else:
        prop_identity = 0
        print("lengths don't match")
    return prop_identity


def parallel_align(input_list=None, n_procs=5):
    """align a certain number of sequences, in parallel, using the multiprocessing module"""
    from multiprocessing import Pool
    import time

    print("DEBUG: we're in parallel_align")
    n_comparisons = len(input_list)
    print("DEBUG: we're about to init time")
    init_time = time.time()
    print("DEBUG: we're about to printout the n of comparisons")
    print(f"Generating {n_comparisons} alignments on {n_procs} processors; this should take ~{round(n_comparisons / n_procs / 100)} seconds.")
    print("DEBUG: we're about to enter the with statement")
    with Pool(n_procs) as p:
        print("DEBUG: we're in the with statement")
        list_results = [i for i in  p.map(align_2_seqs, input_list)]
        print("DEBUG: we're probably done aligning")
    print("DEBUG: we're out of the with statement")
    print(f"Done in {round(time.time() - init_time)} seconds")
    
    return list_results


def is_identical_inclusive(nt1, nt2):
    """
    returns true if the 2 nts are deemed identical:
    are the same nt, or are ambiguous bases that share a possible nt
    """
    dict_IUPAC = {
        "M": ("A", "C"),
        "V": ("A", "C", "G"),
        "H": ("A", "C", "T"),
        "N": ("A", "C", "G", "T"),
        "R": ("A", "G"),
        "D": ("A", "G", "T"),
        "W": ("A", "T"),
        "S": ("C", "G"),
        "B": ("C", "G", "T"),
        "Y": ("C", "T"),
        "K": ("G", "T"),
        "A": ("A",),
        "T": ("T",),
        "G": ("G",),
        "C": ("C",),
        "-": ("-",),
        ".": (".",),
    }
    if any([i in dict_IUPAC[nt2] for i in dict_IUPAC[nt1]]):
        return True
    else:
        return False


def count_identical_bases(target, query):
    """Counts the number of nts identical between 2 aligned sequences
    Indels are considered differences; ambiguous bases are considered identical if one of the possible nts is a match.
    """
    return sum([1 for index in range(len(target)) if is_identical_inclusive(target[index], query[index])])


