import os
from collections import defaultdict
species = ['chicken', 'pig', 'cattle', 'sheep', 'dog']
blast_results = defaultdict(dict)

def parse_blast_results(directory):
    for query_sp in species:
        for ref_sp in species:
            if query_sp != ref_sp: 
                file_name = f"{ref_sp}_ref_{query_sp}_query.txt"
                file_path = os.path.join(directory, file_name)
                if os.path.exists(file_path):
                    with open(file_path, 'r') as f:
                        for line in f:
                            fields = line.strip().split('\t')
                            query_id, ref_id = fields[0], fields[1]
                            blast_results[(query_sp, ref_sp)][query_id] = ref_id
                else:
                    print(f"Warning: File {file_name} not found.")

def find_orthologous_groups():
    orthologs = []
    for chicken_gene in blast_results.get(('chicken', species[1]), {}):
        candidate_group = {'chicken': chicken_gene}
        valid = True
        for i, query_sp in enumerate(species):
            for ref_sp in species[i+1:] + species[:i]:
                if query_sp != ref_sp:
                    query_gene = candidate_group.get(query_sp)
                    if not query_gene:
                        continue
                    if query_gene not in blast_results[(query_sp, ref_sp)]:
                        valid = False
                        break
                    ref_gene = blast_results[(query_sp, ref_sp)][query_gene]
                    if ref_sp not in candidate_group:
                        candidate_group[ref_sp] = ref_gene
                    else:
                        if candidate_group[ref_sp] != ref_gene:
                            valid = False
                            break
                    reverse_query = blast_results.get((ref_sp, query_sp), {}).get(ref_gene)
                    if reverse_query != query_gene:
                        valid = False
                        break
            if not valid:
                break
        if valid and len(candidate_group) == len(species):
            orthologs.append(candidate_group)
    return orthologs

def main(blast_dir, output_file):
    parse_blast_results(blast_dir)
    orthologs = find_orthologous_groups()
    with open(output_file, 'w') as out:
        out.write("chicken\tpig\tcattle\tsheep\tdog\n")
        for group in orthologs:
            out.write(f"{group['chicken']}\t{group['pig']}\t{group['cattle']}\t{group['sheep']}\t{group['dog']}\n")
    print(f"Found {len(orthologs)} orthologous groups. Results written to {output_file}")

blast_directory = "../result" 
output_file = "../result/5species_orthologs.txt"
main(blast_directory, output_file)
