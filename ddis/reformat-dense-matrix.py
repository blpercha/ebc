from collections import defaultdict

def main():
    f = open("../resources/matrix-ebc-paper-dense.tsv", "r")

    gene_connection_map = defaultdict(set)

    for line in f:
        sl = line.split("\t")
        col1_names = sl[0].strip("()").split(",")
        drug = col1_names[0]
        gene = col1_names[1]
        dependency_path = sl[2]
        gene_connection_map[gene].add((drug, dependency_path))
        val = float(sl[4].strip())
    f.close()

    # write new matrix
    f3d = open("/Users/beth/Documents/phd/ebc-ddis/matrix-ebc-paper-ddi-format.tsv", "w")

    for gene in gene_connection_map:
        for drug_path_pair1 in gene_connection_map[gene]:
            drug1 = drug_path_pair1[0]
            path1 = drug_path_pair1[1]
            for drug_path_pair2 in gene_connection_map[gene]:
                drug2 = drug_path_pair2[0]
                path2 = drug_path_pair2[1]
                if drug1 == drug2:
                    continue
                print >> f3d, drug1 + "," + drug2 + "\t" + drug1 + "\t" + drug2 + "\t" + \
                              str(path1) + "\t" + str(path2) + "\t" + gene + "\t1"

    f3d.flush()
    f3d.close()

if __name__ == "__main__":
    main()