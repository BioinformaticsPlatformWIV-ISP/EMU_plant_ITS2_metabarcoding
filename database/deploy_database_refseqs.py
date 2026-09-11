"""
Extracts all fasta records into 1 record per file for ease of use for alignment during pipeline run.
To be replaced with a on-the-fly extraction of the fasta sequences of target species, which should be easy.
"""
from Bio import SeqIO
import pathlib

cur_path = pathlib.Path(__file__).parent.resolve()
print(cur_path)

print("extracting and copying database. this may take a few minutes")
count = 0
for record in SeqIO.parse(cur_path / "species_taxid.fasta", "fasta"):
    with open(cur_path / f"refseqs/{record.id}.fasta", "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
        count+=1
    if count % 10000 == 0:
        print(count)