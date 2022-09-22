import os
import sys


# get cis and trans interactions in separate folder and files
def cis_trans(hic_out, out_dir):
    oc_file = out_dir+"/cis"
    if not os.path.exists(oc_file):
        os.makedirs(oc_file)
    ot_file = out_dir+"/trans"
    if not os.path.exists(ot_file):
        os.makedirs(ot_file)

    with open(hic_out, 'r') as f:
        for i in f.readlines():
            c1 = i.split('\t')[1]
            p1 = i.split('\t')[2]
            s1 = i.split('\t')[3]
            c2 = i.split('\t')[4]
            p2 = i.split('\t')[5]
            s2 = i.split('\t')[6]
            if c1 == c2:
                file_e = open(oc_file + "/chr_" + c1 + ".txt", "a")
                file_e.write(c1 + "\t" + p1 + "\t" + s1 + "\t" + c2 + "\t" + p2 + "\t" + s2 + "\n")
                file_e.close()

            if c1 != c2:
                file_e = open(ot_file + "/chr_" + c1 + "_" + c2 + ".txt", "a")
                file_e.write(c1 + "\t" + p1 + "\t" + s1 + "\t" + c2 + "\t" + p2 + "\t" + s2 + "\n")
                file_e.close()


def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)


# get gene name by coordinates from cis and trans files
def get_gene_name_by_cord(gtf_file, chr_name, location, strand):
    gene_name = ""
    chr_name = str(chr_name.strip())
    location = int(location.strip())
    strand = str(strand.strip())
    
    with open(gtf_file) as f:
        for line_s in f.readlines():
            chr_col = str(line_s.split("\t")[0])
            s_gene = int(line_s.split("\t")[3])
            e_gene = int(line_s.split("\t")[4])
            m_strand = str(line_s.split("\t")[6])

            if chr_name == chr_col and s_gene <= location <= e_gene and strand == m_strand:
                gene_name = line_s.split("\t")[8].split(":")[1].split(";")[0]
    return chr_name+"\t"+str(location)+"\t"+strand+"\t"+gene_name


if __name__ == '__main__':
    in_file = sys.argv[1]
    output_dir = sys.argv[2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    cis_trans(in_file, output_dir)

    cis_lst = []
    for files in os.listdir(output_dir+"/cis"):
        files_name = output_dir+"/cis/"+files
        with open(files_name, 'rb') as fp:
            c_generator = _count_generator(fp.raw.read)
            list_count = sum(buffer.count(b'\n') for buffer in c_generator)
            cis_lst.append(list_count)
            cis_count = sum(cis_lst)

    trans_lst = []
    for files in os.listdir(output_dir+"/trans"):
        files_name = output_dir+"/trans/"+files
        with open(files_name, 'rb') as fp:
            c_generator = _count_generator(fp.raw.read)
            list_count = sum(buffer.count(b'\n') for buffer in c_generator)
            trans_lst.append(list_count)
            trans_count = sum(trans_lst)

    ratio = round((cis_count/trans_count)*100, 2)
    print("cis and trans interactions in "+in_file+" is: "+str(cis_count)+" and "+str(trans_count))
    print("cis and trans interaction ratio(%) in "+in_file+" is: "+str(ratio))

    gtf = "/home/ajeet/Projects/PFusion/Arabidopsis_thaliana.TAIR10.51.gene.gff3"
    for files in os.listdir(output_dir+"/cis"):
        files_name = output_dir+"/cis/"+files
        file = open(files_name+".gene.tab", "a")
        with open(files_name, 'r') as fp:
            for lines in fp.readlines():
                chr_1n = lines.split("\t")[0]
                loc_1n = lines.split("\t")[1]
                s_1n = lines.split("\t")[2]
                chr_2n = lines.split("\t")[3]
                loc_2n = lines.split("\t")[4]
                s_2n = lines.split("\t")[5]

                part_a = get_gene_name_by_cord(gtf, chr_1n, loc_1n, s_1n)
                part_b = get_gene_name_by_cord(gtf, chr_2n, loc_2n, s_2n)
                # print(part_a.strip(), part_b.strip(), sep='\t')

                if part_a and part_b is not None:
                    file.write(part_a.strip()+"\t"+part_b.strip()+"\n")
        file.close()

    for files in os.listdir(output_dir+"/trans"):
        files_name = output_dir+"/trans/"+files
        file = open(files_name+".gene.tab", "a")
        with open(files_name, 'r') as fp:
            for lines in fp.readlines():
                chr_1n = lines.split("\t")[0]
                loc_1n = lines.split("\t")[1]
                s_1n = lines.split("\t")[2]
                chr_2n = lines.split("\t")[3]
                loc_2n = lines.split("\t")[4]
                s_2n = lines.split("\t")[5]

                part_a = get_gene_name_by_cord(gtf, chr_1n, loc_1n, s_1n)
                part_b = get_gene_name_by_cord(gtf, chr_2n, loc_2n, s_2n)
                # print(part_a.strip(), part_b.strip(), sep='\t')

                if part_a and part_b is not None:
                    file.write(part_a.strip()+"\t"+part_b.strip()+"\n")
        file.close()
