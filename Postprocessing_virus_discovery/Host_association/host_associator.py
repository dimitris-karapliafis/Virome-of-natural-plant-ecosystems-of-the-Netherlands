
import pandas as pd
from sys import argv
from ete3 import NCBITaxa
import numpy as np

def read_annot_file(annot_file: str) -> pd.DataFrame:
    """
    Read the annotation file and return a pandas DataFrame
    :param annot_file:
    :return:
    """
    return pd.read_csv(annot_file, sep='\t',header=None)


def read_mmseqs_tax_file(h_assoc_file: str) -> pd.DataFrame:
    """
    Read the h_assoc file and return a pandas DataFrame
    :param h_assoc_file:
    :return:
    """
    df = pd.read_csv(h_assoc_file, sep='\t')
    df = df.replace("#NAME?", "NA")

    return df


def file_tester(file):

    with open(file) as fn:
        for line in fn:
            line.count("\t")
            if line.count("\t") != 18:
                print(line.split("\t"))
            else:
                continue

def subset_df(df, cols):
    """
    Subset a DataFrame by selecting only specific columns based on their index
    :param df:
    :param cols: list of column indices
    """
    return df.iloc[:, cols]



def merge_dfs(df1, df2, on):
    """
    Merge two DataFrames
    :param df1:
    :param df2:
    :return:
    """
    return pd.merge(df1, df2, how='outer', on= on)


def filter_nones(x):
    return [i for i in x if i != None]

def get_virus_family_genus_species(virus_tax_id):
    ncbi = NCBITaxa()
    try:
        lineage = ncbi.get_lineage(virus_tax_id)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        family = 'NA'
        genus = 'NA'
        species = 'NA'
        for taxid in lineage:
            if ranks[taxid] == 'family':
                family = names[taxid]
            if ranks[taxid] == 'genus':
                genus = names[taxid]
            if ranks[taxid] == 'species':
                species = names[taxid]
        return family, genus, species

    except ValueError:
        return 'NA', 'NA', 'NA'



def parse_mmseqs_taxonomy(mmseqs_taxonomy_path, family_level_dict, genus_level_dict, species_level_dict, out_fn):


    with open(mmseqs_taxonomy_path, 'r') as f, open(out_fn, 'w') as out_f:
        out_f.write("vOTU\ttaxonomy_id\ttaxonomic_rank\ttaxonomy_name\tNo. of fragments retained"
                    "\tNo. of fragments taxonomically assigned\tNo. of fragments in agreement with taxonomic label"
                    "\tSupport_received\tLineage\tFamily_level_host_association\tGenus_level_host_association"
                    "\tSpecies_level_host_association\n")

        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            virus_tax_id = line[1]
            family, genus, species = get_virus_family_genus_species(virus_tax_id)

            if family in family_level_dict:
                no_none_family = filter_nones(family_level_dict[family])
                no_none_family = set(no_none_family)
                family_str = ",".join(no_none_family)
                if family == 'NA' or family_str == '':
                    family_str = 'NA'
                    if len(line) == 9:
                        line.append(family_str)
                    else:
                        line.append("NA")
                        line.append(family_str)
                else:
                    if len(line) == 9:
                        line.append(family_str)
                    else:
                        line.append("NA")
                        line.append(family_str)
            else:
                if len(line) == 9:
                    line.append("NA")
                else:
                    line.append("NA")
                    line.append('NA')

            if genus in genus_level_dict:
                no_none_genus = filter_nones(genus_level_dict[genus])
                no_none_genus = set(no_none_genus)
                genus_str = ",".join(no_none_genus)
                if genus == 'NA' or genus_str == '':
                    genus_str = 'NA'
                    if len(line) == 10:
                        line.append(genus_str)
                    else:
                        line.append("NA")
                        line.append(genus_str)
                else:
                    if len(line) == 10:
                        line.append(genus_str)
                    else:
                        line.append("NA")
                        line.append(genus_str)
            else:
                if len(line) == 10:
                    line.append('NA')
                else:
                    line.append("NA")
                    line.append('NA')

            if species in species_level_dict:

                no_none_species = filter_nones(species_level_dict[species])
                no_none_species = set(no_none_species)
                species_str = ",".join(no_none_species)
                if species == 'NA' or species_str == '':
                    species_str = 'NA'
                    if len(line) == 11:
                        line.append(species_str)
                    else:
                        line.append("NA")
                        line.append(species_str)
                else:
                    if len(line) == 11:
                        line.append(species_str)
                    else:
                        line.append("NA")
                        line.append(species_str)

            else:
                if len(line) == 11:
                    line.append('NA')
                else:
                    line.append("NA")
                    line.append('NA')

            out_f.write("\t".join(line) + "\n")

    return out_fn


def ictv_mmseqs_tax_wrapper(df, mmseqs_taxonomy_path, out_fn):
    species_level_dict = {}
    genus_level_dict = {}
    family_level_dict = {}

    for index, row in df.iterrows():
        species = row['Species']
        genus = row['Genus']
        family = row['Family']
        source = row['Host source']
        for s in source.split(","):
            new_source = s.strip()


            if species not in species_level_dict:
                species_level_dict[species] = []
            species_level_dict[species].append(new_source)

            if genus not in genus_level_dict:
                genus_level_dict[genus] = []
            genus_level_dict[genus].append(new_source)

            if family not in family_level_dict:
                family_level_dict[family] = []
            family_level_dict[family].append(new_source)

    out_fn = parse_mmseqs_taxonomy(mmseqs_taxonomy_path, family_level_dict, genus_level_dict, species_level_dict,
                                   out_fn)

    return out_fn


def read_xlsx_file(xlsx_file):
    df = pd.read_excel(xlsx_file)
    return df


def extract_tax_levels(df):

    spec_set = set()
    gen_set = set()
    fam_set = set()
    ord_set = set()
    cls_set = set()
    phy_set = set()
    kin_set = set()
    sup_set = set()

    tmp_spec = df['Species']
    tmp_gen = df['Genus']
    tmp_fam = df['Family']
    tmp_ord = df['Order']
    tmp_cls = df['Class']
    tmp_phy = df['Phylum']
    tmp_kin = df['Kingdom']



    for spec, gen, fam, ord_, cls, phy, kin  in zip(tmp_spec,
                                                         tmp_gen,
                                                         tmp_fam,
                                                         tmp_ord,
                                                         tmp_cls,
                                                         tmp_phy,
                                                         tmp_kin):

        spec_set.add(spec)
        gen_set.add(gen)
        fam_set.add(fam)
        ord_set.add(ord_)
        cls_set.add(cls)
        phy_set.add(phy)
        kin_set.add(kin)

    return spec_set, gen_set, fam_set, ord_set, cls_set, phy_set, kin_set


def taxonomizer_2(merged_df, ictv_xlsx_file):

    taxonomizer_dict = {}

    species_set, genus_set, family_set, order_set, class_set, phylum_set, kingdom_set = (
        extract_tax_levels(read_xlsx_file(ictv_xlsx_file)))

    for _,row in merged_df.iterrows():
        votu = row.vOTU
        if type(row.genomad_taxonomy) == str:
            genomad_taxonomy = row.genomad_taxonomy.split(";")
        else:
            genomad_taxonomy = ["NA"]
        if type(row.palmdb_taxonomy) == str:
            palmdb_taxonomy = row.palmdb_taxonomy.split(";")
        else:
            palmdb_taxonomy = ["NA"]
        geno_tax = {}
        palmdb_tax = {}
        plant_vir_db_tax = {}
        mmseqs_tax = {}
        for element in ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]:
            geno_tax[element] = "NA"
            palmdb_tax[element] = "NA"
            plant_vir_db_tax[element] = "NA"
            mmseqs_tax[element] = "NA"

        for tax in genomad_taxonomy:

            if tax in kingdom_set:
                geno_tax["kingdom"] = tax
            if tax in phylum_set:
                geno_tax["phylum"] = tax
            if tax in class_set:
                geno_tax["class"] = tax
            if tax in order_set:
                geno_tax["order"] = tax
            if tax in family_set:
                geno_tax["family"] = tax
            if tax in genus_set:
                geno_tax["genus"] = tax
            if tax in species_set:
                geno_tax["species"] = tax



        for tax in palmdb_taxonomy:

            if tax in kingdom_set:
                palmdb_tax["kingdom"] = tax
            if tax in phylum_set:
                palmdb_tax["phylum"] = tax
            if tax in class_set:
                palmdb_tax["class"] = tax
            if tax in order_set:
                palmdb_tax["order"] = tax
            if tax in family_set:
                palmdb_tax["family"] = tax
            if tax in genus_set:
                palmdb_tax["genus"] = tax
            if tax in species_set:
                palmdb_tax["species"] = tax

        if type(row.plant_virus_db_taxonomy) == str:
            plant_vir_db_taxonomy = row.plant_virus_db_taxonomy.split("_")
        else:
            plant_vir_db_taxonomy = ["NA"]
        if "NA" not in plant_vir_db_taxonomy:
            plant_vir_db_tax["family"] = plant_vir_db_taxonomy[0]

        mmseqs_k = row.kingdom
        mmseqs_p = row.phylum
        mmseqs_c = row["class"]
        mmseqs_o = row.order
        mmseqs_f = row.family
        mmseqs_g = row.genus
        mmseqs_s = row.species
        if type(mmseqs_k) == str and mmseqs_k != "Not annotated by mmseqs easy-taxonomy":
            mmseqs_tax["kingdom"] = mmseqs_k
        else:
            mmseqs_tax["kingdom"] = "NA"
        if type(mmseqs_p) == str and mmseqs_p != "Not annotated by mmseqs easy-taxonomy":
            mmseqs_tax["phylum"] = mmseqs_p
        else:
            mmseqs_tax["phylum"] = "NA"
        if type(mmseqs_c) == str and mmseqs_c != "Not annotated by mmseqs easy-taxonomy":
            mmseqs_tax["class"] = mmseqs_c
        else:
            mmseqs_tax["class"] = "NA"
        if type(mmseqs_o) == str and mmseqs_o != "Not annotated by mmseqs easy-taxonomy":
            mmseqs_tax["order"] = mmseqs_o
        else:
            mmseqs_tax["order"] = "NA"
        if type(mmseqs_f) == str and mmseqs_f != "Not annotated by mmseqs easy-taxonomy":
            mmseqs_tax["family"] = mmseqs_f
        else:
            mmseqs_tax["family"] = "NA"
        if type(mmseqs_g) == str and mmseqs_g != "Not annotated by mmseqs easy-taxonomy":
            mmseqs_tax["genus"] = mmseqs_g
        else:
            mmseqs_tax["genus"] = "NA"
        if type(mmseqs_s) == str and mmseqs_s != "Not annotated by mmseqs easy-taxonomy":
            mmseqs_tax["species"] = mmseqs_s
        else:
            mmseqs_tax["species"] = "NA"

        taxonomizer_dict[votu] = {"genomad": geno_tax, "palmdb": palmdb_tax,
                                  "plant_virus_db": plant_vir_db_tax, "mmseqs": mmseqs_tax}
    return taxonomizer_dict


def parse_ictv_xlsx_file(df):
    species_level_dict = {}
    genus_level_dict = {}
    family_level_dict = {}

    for index, row in df.iterrows():
        species = row['Species']
        genus = row['Genus']
        family = row['Family']
        source = row['Host source']
        for s in source.split(","):
            new_source = s.strip()


            if species not in species_level_dict:
                species_level_dict[species] = set()
            species_level_dict[species].add(new_source)

            if genus not in genus_level_dict:
                genus_level_dict[genus] = set()
            genus_level_dict[genus].add(new_source)

            if family not in family_level_dict:
                family_level_dict[family] = set()
            family_level_dict[family].add(new_source)

    return species_level_dict, genus_level_dict, family_level_dict



def find_deeper_annotation(merged_df, ictv_xlsx_file):

    tax_df = pd.DataFrame()
    taxonomizer_dict = taxonomizer_2(merged_df, ictv_xlsx_file)
    votus = merged_df.vOTU
    d_annot_dict = {}
    consensus_species = []
    consensus_genus = []
    consensus_family = []
    consensus_order = []
    consensus_class = []
    consensus_phylum = []
    consensus_kingdom = []
    consensus_superkingdom = []
    species_lvl_host = []
    genus_lvl_host = []
    family_lvl_host = []
    best_host_as = []
    best_host_as_lvl = []
    host_as_tax_group = []

    tool_list = []

    species_level_dict, genus_level_dict, family_level_dict = parse_ictv_xlsx_file(read_xlsx_file(ictv_xlsx_file))


    for _,row in merged_df.iterrows():

        votu = row.vOTU
        species_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['species']}",
                       f"palmdb: {taxonomizer_dict[votu]['palmdb']['species']}",
                       f"plant_virus_db: {taxonomizer_dict[votu]['plant_virus_db']['species']}",
                       f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['species']}"]
        genus_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['genus']}",
                     f"palmdb: {taxonomizer_dict[votu]['palmdb']['genus']}",
                     f"plant_virus_db: {taxonomizer_dict[votu]['plant_virus_db']['genus']}",
                     f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['genus']}"]
        family_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['family']}",
                      f"palmdb: {taxonomizer_dict[votu]['palmdb']['family']}",
                      f"plant_virus_db: {taxonomizer_dict[votu]['plant_virus_db']['family']}",
                      f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['family']}"]
        order_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['order']}",
                     f"palmdb: {taxonomizer_dict[votu]['palmdb']['order']}",
                     f"plant_virus_db: NA", f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['order']}"]
        class_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['class']}",
                     f"palmdb: {taxonomizer_dict[votu]['palmdb']['class']}",
                     f"plant_virus_db: NA", f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['class']}"]
        phylum_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['phylum']}",
                      f"palmdb: {taxonomizer_dict[votu]['palmdb']['phylum']}",
                      f"plant_virus_db: NA", f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['phylum']}"]
        kingdom_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['kingdom']}",
                       f"palmdb: {taxonomizer_dict[votu]['palmdb']['kingdom']}",
                       f"plant_virus_db: NA", f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['kingdom']}"]

        superkingdom_lvl = [f"genomad: {taxonomizer_dict[votu]['genomad']['superkingdom']}",
                            f"palmdb: {taxonomizer_dict[votu]['palmdb']['superkingdom']}",
                            f"plant_virus_db: NA", f"mmseqs: {taxonomizer_dict[votu]['mmseqs']['superkingdom']}"]

        consensus_species.append(f"{', '.join(species_lvl)}")
        consensus_genus.append(f"{', '.join(genus_lvl)}")
        consensus_family.append(f"{', '.join(family_lvl)}")
        consensus_order.append(f"{', '.join(order_lvl)}")
        consensus_class.append(f"{', '.join(class_lvl)}")
        consensus_phylum.append(f"{', '.join(phylum_lvl)}")
        consensus_kingdom.append(f"{', '.join(kingdom_lvl)}")
        consensus_superkingdom.append(f"{', '.join(superkingdom_lvl)}")

        species_lvl_host.append(", ".join(ictv_host_checker(species_level_dict, species_lvl)))
        genus_lvl_host.append(",".join(ictv_host_checker(genus_level_dict, genus_lvl)))
        family_lvl_host.append(", ".join(ictv_host_checker(family_level_dict, family_lvl)))
        bh, bh_lvl = best_host_picker_2(ictv_host_checker(family_level_dict, family_lvl),
                                        ictv_host_checker(genus_level_dict,genus_lvl),
                                        ictv_host_checker(species_level_dict, species_lvl))

        best_host_as.append("|".join(bh))
        best_host_as_lvl.append(bh_lvl)
        tax_group = set()
        tool_set = set()

        if bh_lvl == "species":
            for element in species_lvl:
                tool, tax = element.split(":")
                if tax.strip() != "NA" and type(tax) == str:
                    tax_group.add(tax.strip())
                    tool_set.add(tool.strip())
        elif bh_lvl == "genus":
            for element in genus_lvl:
                tool, tax = element.split(":")
                if tax.strip() != "NA" and type(tax) == str:
                    tax_group.add(tax.strip())
                    tool_set.add(tool.strip())
        elif bh_lvl == "family":
            for element in family_lvl:
                tool, tax = element.split(":")
                if tax.strip() != "NA" and type(tax) == str:
                    tax_group.add(tax.strip())
                    tool_set.add(tool.strip())
        else:
            tax_group.add("NA")
            tool_set.add("NA")
        host_as_tax_group.append(", ".join(tax_group))
        tool_list.append(", ".join(tool_set))

    best_host_as_host = []
    for element in best_host_as:
        if len(element.split("|")) > 1:
            best_host_as_host.append("Conflict! Check manually")
        else:
            if element != "NA":
                best_host_as_host.append(element.split(":")[-1])
            else:
                best_host_as_host.append("NA")


    tax_df = pd.DataFrame({"vOTU": votus,"combined_family_annotation": consensus_family,
                           "combined_genus_annotation": consensus_genus, "combined_species_annotation": consensus_species,
                           "ictv_host_family (combined)": family_lvl_host, "ictv_host_genus (combined)": genus_lvl_host,
                           "ictv_host_species (combined)": species_lvl_host,
                           "best_host_association_overview (tool:host range)": best_host_as,
                           "approach_of_taxonomic_annotation": tool_list,
                           "taxonomic_level_for_host_association": best_host_as_lvl,
                           "taxonomic_group_for_host_association": host_as_tax_group,
                           "host_association": best_host_as_host
                          })


    tax_df.to_csv("host_association_interm_ictv.tsv", sep="\t", index=False)

    return tax_df



def ictv_host_checker(tax_lvl_dict, tax_lvl_list):

    host_list = []
    for element in tax_lvl_list:
        tool, tax = element.split(":")
        tax = tax.strip()
        if tax in tax_lvl_dict:

            host_list.append(f"{tool}:{', '.join(tax_lvl_dict[tax])}" )
        else:
            host_list.append(f"{tool}: NA")
    return host_list


def best_host_picker_2(family_lvl_host, genus_lvl_host, species_lvl_host):
    interm_host_list = []
    best_host_set = set()
    species = False
    genus = False
    family = False

    for element in species_lvl_host:
        tool, host = element.split(":")
        host = host.strip()
        if host != "NA":
            species = True
            interm_host_list.append(element)
    if len(interm_host_list) > 1:
        tools = []
        hosts = []
        for element in interm_host_list:
            tool, host = element.split(":")
            tools.append(tool)
            hosts.append(host)
        # If the host annotations are the same for any tools, then create a concatanated string
        # with the tools with the same annotation and the annotation

        for host in set(hosts):
            if len(set(hosts)) > 1:
                best_host_set.add(f"{','.join([tools[i] for i, x in enumerate(hosts) if x == host])}:{host}")
            else:
                best_host_set.add(f"{tools[hosts.index(host)]}:{host}")
    elif len(interm_host_list) == 1:
        best_host_set.add(interm_host_list[0])

    if species!= True and not best_host_set:
        for element in genus_lvl_host:
            tool, host = element.split(":")
            host = host.strip()
            if host != "NA":
                genus = True
                interm_host_list.append(element)
        if len(interm_host_list) > 1:
            tools = []
            hosts = []
            for element in interm_host_list:
                tool, host = element.split(":")
                tools.append(tool)
                hosts.append(host)

            # If the host annotations are the same for any tools, then create a concatanated string
            # with the tools with the same annotation and the annotation
            for host in set(hosts):
                if len(set(hosts)) > 1:
                    best_host_set.add(f"{','.join([tools[i] for i, x in enumerate(hosts) if x == host])}:{host}")
                else:
                    best_host_set.add(f"{tools[hosts.index(host)]}:{host}")

        elif len(interm_host_list) == 1:
            best_host_set.add(interm_host_list[0])


    if genus!= True and not best_host_set:
        for element in family_lvl_host:
            tool, host = element.split(":")
            host = host.strip()
            if host != "NA":
                family = True
                interm_host_list.append(element)
        if len(interm_host_list) > 1:
            tools = []
            hosts = []
            for element in interm_host_list:
                tool, host = element.split(":")
                tools.append(tool)
                hosts.append(host)
            # If the host annotations are the same for any tools, then create a concatanated string
            # with the tools with the same annotation and the annotation
            for host in set(hosts):
                if len(set(hosts)) > 1:
                    best_host_set.add(f"{','.join([tools[i] for i, x in enumerate(hosts) if x == host])}:{host}")
                else:
                    best_host_set.add(f"{tools[hosts.index(host)]}:{host}")
        elif len(interm_host_list) == 1:
            best_host_set.add(interm_host_list[0])


    if not species and not genus and not family:
        best_host_set.add("NA")
        return list(best_host_set), "NA"
    elif species:
        return list(best_host_set), "species"
    elif genus:
        return list(best_host_set), "genus"
    elif family:
        return list(best_host_set), "family"


def get_virus_lineage(virus_tax_id):
    ncbi = NCBITaxa()
    try:
        lineage = ncbi.get_lineage(virus_tax_id)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        superkingdom = 'NA'
        kingdom = 'NA'
        phylum = 'NA'
        class_= 'NA'
        order = 'NA'
        family = 'NA'
        genus = 'NA'
        species = 'NA'
        for taxid in lineage:
            if ranks[taxid] == 'superkingdom':
                superkingdom = names[taxid]
            if ranks[taxid] == 'kingdom':
                kingdom = names[taxid]
            if ranks[taxid] == 'phylum':
                phylum = names[taxid]
            if ranks[taxid] == 'class':
                class_ = names[taxid]
            if ranks[taxid] == 'order':
                order = names[taxid]
            if ranks[taxid] == 'family':
                family = names[taxid]
            if ranks[taxid] == 'genus':
                genus = names[taxid]
            if ranks[taxid] == 'species':
                species = names[taxid]
        return superkingdom, kingdom, phylum, class_, order, family, genus, species

    except ValueError:
        return "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"


def write_mmseqs_file(mmseqs_ictv_tax_fn, annotated_fn, out_fn):
    mmseqs_votus = []
    with open(mmseqs_ictv_tax_fn) as in_handle, open(out_fn, "w") as out_handle, open(annotated_fn) as annot_handle:
        for line in in_handle:
            if line.startswith("vOTU"):
                out_handle.write("vOTU\t"
                                 "vOTU_length\t"
                                 "virus_tax_id\t"
                                 "taxonomic_rank\t"
                                 "taxonomy_name\t"
                                 "superkingdom\t"
                                 "kingdom\t"
                                 "phylum\t"
                                 "class\t"
                                 "order\t"
                                 "family\t"
                                 "genus\t"
                                 "species\t"
                                 "family_level_host_ictv\t"
                                 "genus_level_host_ictv\t"
                                 "species_level_host_ictv\n")

                continue
            line = line.strip().split("\t")
            votu = line[0]
            votu_len = votu.split("_")[-3]
            votu_taxid = line[1]
            tax_rank = line[2]
            tax_name = line[3]
            family_level_host_ictv = line[9]
            genus_level_host_ictv = line[10]
            species_level_host_ictv = line[11]
            domain, kingdom, phylum, class_, order, family, genus, species = get_virus_lineage(votu_taxid)
            mmseqs_votus.append(votu)
            out_handle.write(votu + "\t" +
                             votu_len + "\t"
                             + votu_taxid + "\t"
                             + tax_rank + "\t"
                             + tax_name + "\t"
                             + domain + "\t"
                             + kingdom + "\t"
                             + phylum + "\t"
                             + class_ + "\t"
                             + order + "\t"
                             + family + "\t"
                             + genus + "\t"
                             + species + "\t"
                             + family_level_host_ictv + "\t"
                             + genus_level_host_ictv + "\t"
                             + species_level_host_ictv + "\n")
        for annot_line in annot_handle:
            annot_line = annot_line.strip().split("\t")
            if annot_line[0] not in mmseqs_votus:
                non_mmseqs_votu_str = "Not annotated by mmseqs easy-taxonomy"
                votu_length = annot_line[0].split("_")[-3]
                out_handle.write(annot_line[0] + "\t" + votu_length + "\t" + (non_mmseqs_votu_str + "\t") * 13 + "\n")
    return out_fn


if __name__ == "__main__":
    annot_file = argv[1]
    mmseqs_tax_file = argv[2]
    ictv_xlsx_file = argv[3]
    annot_df = read_annot_file(annot_file)
    annot_subset = annot_df.iloc[:, [0, 5, 7, 8]]
    title_line_annot = ["vOTU", "genomad_taxonomy", "palmdb_taxonomy", "plant_virus_db_taxonomy"]
    annot_subset.columns = title_line_annot
    mmseqs_df = read_mmseqs_tax_file(mmseqs_tax_file)
    title_line_mmseqs = ["vOTU","mmseqs_easy_taxonomy_lineage"]
    mmseqs_subset = mmseqs_df.iloc[:, [0, 8]]
    mmseqs_subset.columns = title_line_mmseqs
    ictv_df = read_xlsx_file(ictv_xlsx_file)
    mmseqs_ictv_tax_fn = ictv_mmseqs_tax_wrapper(ictv_df, mmseqs_tax_file, "mmseqs_ictv_merged.tsv")
    h_assoc_file = write_mmseqs_file(mmseqs_ictv_tax_fn, annot_file, "mmseqs_ictv_host_association.tsv")
    h_assoc_df = read_mmseqs_tax_file(h_assoc_file)

    merged_df = merge_dfs(annot_subset, h_assoc_df, on="vOTU")
    tax_df = find_deeper_annotation(merged_df, ictv_xlsx_file)
    merged_df = merge_dfs(annot_subset, mmseqs_subset, on="vOTU")
    tax_subset = subset_df(tax_df, [0,7, 8, 9, 10, 11])
    merged_df = merge_dfs(merged_df, tax_subset, on="vOTU")
    merged_df.fillna("NA", inplace=True)
    merged_df.to_csv("host_association_ictv.tsv", sep="\t", index=False)
    plant_subset_df = merged_df[merged_df['host_association'] == "plants"]
    plant_subset_df.to_csv("plant_host_association_ictv.tsv", sep="\t", index=False)

