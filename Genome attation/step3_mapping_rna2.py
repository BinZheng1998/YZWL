# replace_protein_id_with_gene_id_in_gff.py

import sys

def load_mapping(mapping_file):
    """
    Load mapping file into a dictionary where protein ID maps to gene ID.
    """
    protein_to_gene = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            if line.strip():  # skip empty lines
                parts = line.strip().split('\t')
                protein_id = parts[0]
                gene_id = parts[1]
                protein_to_gene[protein_id] = gene_id
    return protein_to_gene

def process_gff(mapping_file, gff_file):
    """
    Process GFF file and replace protein IDs with 'rna-gene_id' in 'transdecoder' lines.
    """
    protein_to_gene = load_mapping(mapping_file)
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # skip comment lines
                print(line.strip())
                continue
            
            parts = line.strip().split('\t')
            
            if len(parts) < 9:
                # If there are less than 9 columns, skip this line
                print(line.strip())
                continue
            
            if parts[1] == 'transdecoder':
                # Process 'transdecoder' lines
                attributes = parts[8].split(';')
                new_attributes = []
                
                for attr in attributes:
                    if attr.startswith('ID='):
                        protein_id = attr.split('=')[1]
                        if protein_id in protein_to_gene:
                            gene_id = protein_to_gene[protein_id]
                            new_attributes.append(f"ID=rna-predicted-{gene_id}")
                        else:
                            new_attributes.append(attr)
                    elif attr.startswith('Parent='):
                        protein_id = attr.split('=')[1]
                        if protein_id in protein_to_gene:
                            gene_id = protein_to_gene[protein_id]
                            new_attributes.append(f"Parent=rna-predicted-{gene_id}")
                        else:
                            new_attributes.append(attr)
                    else:
                        new_attributes.append(attr)
                
                parts[8] = ';'.join(new_attributes)
            
            # Output the modified line
            print('\t'.join(parts))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python replace_protein_id_with_gene_id_in_gff.py mapping_file gff_file > output_file")
        sys.exit(1)
    
    mapping_file = sys.argv[1]
    gff_file = sys.argv[2]
    
    process_gff(mapping_file, gff_file)
