import collections

from Bio.KEGG import REST

import networks


def query_KEGG(database="pathway", organism="sce"):
    raw_module_data = REST.kegg_list(database, organism).read()

    module_data = []
    for line in raw_module_data.rstrip().split("\n"):
        entry, description = line.split("\t")
        module_data.append((entry, description.replace('/', '|')))

    # Get the genes for modules and add them to a list
    module_genes_dict = collections.defaultdict(list)
    for module in module_data:
        entry, description = module
        module_file = REST.kegg_get(entry).read()  # query and read each module

        # iterate through each KEGG file, keeping track of which section
        # of the file we're in, only read the gene in each module
        current_section = None
        for line in module_file.rstrip().split("\n"):
            section = line[:12].strip()  # section names are within 12 columns
            if not section == "":
                current_section = section
            if current_section == "GENE":
                gene_name = line[12:].split(' ')[0]
                module_genes_dict[description].append(networks.convert_to_std(gene_name))

    for description in module_genes_dict.keys():
        module_genes_dict[description] = list(set(module_genes_dict[description]))
    return module_genes_dict
