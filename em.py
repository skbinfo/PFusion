#! /usr/bin/env python3
import sys
fuma_file = sys.argv[1]
gtf_file = sys.argv[2]
output_file = sys.argv[3]

def get_gene(gtf, g, chr, loc):
    a = ""
    with open(gtf) as f:
        for line_s in f.readlines():
            if not line_s.startswith("#"):
                if g in line_s and line_s.split("\t")[2] == "exon":
                    chr_col = str(line_s.split("\t")[0].strip())
                    s_gene = int(line_s.split("\t")[3].strip())
                    e_gene = int(line_s.split("\t")[4].strip())
                    if chr == chr_col:
                        if loc-1 == s_gene or loc == s_gene or loc+1 == s_gene or loc-1 == e_gene or loc == e_gene or loc+1 == e_gene:
                            a = g+"\t"+chr+"\t"+str(loc)+"\tE"
                            break
                        elif s_gene+1 < loc < e_gene-1:
                            a = g+"\t"+chr+"\t"+str(loc)+"\tM"
                            break
                        else:
                            a = g+"\t"+chr+"\t"+str(loc)+"\tU"
    return a

if __name__ == '__main__':
    with open(fuma_file, 'r') as fp:
        for lines in fp.readlines():
            g1_g2  = lines.split("\t")[0]
            chr_1n = str(lines.split("\t")[1].strip())
            loc_1n = int(lines.split("\t")[2].strip())
            chr_2n = str(lines.split("\t")[3].strip())
            loc_2n = int(lines.split("\t")[4].strip())

            if g1_g2.count("gene") == 2:
                g1 = g1_g2.split("_")[0].split(":")[1]
                g2 = g1_g2.split("_")[1].split(":")[1]
                
                em1 = get_gene(gtf_file, g1, chr_1n, loc_1n)
                em2 = get_gene(gtf_file, g2, chr_2n, loc_2n)

                if em1 and em2 is not None:
                    file = open(output_file, "a")
                    file.write(g1_g2+"\t"+em1+"\n")
                    file.write(g1_g2+"\t"+em2+"\n")
                    file.close()
                elif em1 is None and em2 is not None:
                    print(g1_g2, chr_1n, loc_1n, " check!! ")
                elif em2 is None and em1 is not None:
                    print(g1_g2, chr_2n, loc_2n, " check!! ")
                else:
                    print(g1_g2, chr_1n, loc_1n, chr_2n, loc_2n, " check!! ")
                
            elif g1_g2.count("gene") > 2:
                all_gene1 = g1_g2.split("_")[0]
                all_gene2 = g1_g2.split("_")[1]
                if all_gene1.count("gene") > 1:
                    anum = all_gene1.count("gene")
                    for num in range(1, anum + anum):
                        if num % 2 != 0:
                            g1a = all_gene1.split(":")[num]
                            em1 = get_gene(gtf_file, g1a, chr_1n, loc_1n)
                            if em1 is not None:
                                file = open(output_file, "a")
                                file.write(g1_g2+"\t"+em1+"\n")
                                file.close()
                            elif em1 is None:
                                print(g1_g2, chr_1n, loc_1n, " check!! ")

                if all_gene2.count("gene") > 1:
                    anum = all_gene2.count("gene")
                    for num in range(1, anum + anum):
                        if num % 2 != 0:
                            g2a = all_gene2.split(":")[num]
                            em2 = get_gene(gtf_file, g2a, chr_2n, loc_2n)
                            if em2 is not None:
                                file = open(output_file, "a")
                                file.write(g1_g2+"\t"+em2+"\n")
                                file.close()
                            elif em2 is None:
                                print(g1_g2, chr_2n, loc_2n, " check!! ")
