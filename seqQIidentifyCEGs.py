'''
@author: abelsousa
'''

import sys, getopt, time

class MyCegs:
    
    def __init__(self, strategy, hmmsearch_output_ref_as, expression_iso, expression_gene, ceg_min_cov, prot_min_cov):
        self.hmmsearch_output_ref_as = self.readHmmsearch(hmmsearch_output_ref_as)
        self.ceg_min_cov = ceg_min_cov
        self.prot_min_cov = prot_min_cov
        if strategy == "reference-based":
            self.isoforms = self.readCuffIso(expression_iso)
            self.genes = self.readCuffGene(expression_gene)
        elif strategy == "denovo":
            self.isoforms = self.readRSEMIso(expression_iso)
            self.genes = self.readRSEMGene(expression_gene)
        else:
            print "This strategy type is not allowed! Please input 'reference-based' or 'denovo' (without quotes)"
            sys.exit()
    


    def readHmmsearch(self, output):
        f = open(output, "r")
        lines = f.readlines()
        f.close()
        line = 0
        query_id = ""
        target_id = ""
        d = {}
        alin = 0
        while line < len(lines):
            if lines[line][0][0] != "#":
                new_line = filter(None, lines[line].strip().split(" "))
                if new_line[3] == query_id:
                    if d.has_key(query_id) and new_line[0] == target_id:
                        d[query_id][alin].append((new_line[0], new_line[2], new_line[5], new_line[12], (new_line[15], new_line[16]), (new_line[17], new_line[18])))
                    elif d.has_key(query_id) and new_line[0] != target_id:
                        d[query_id].append([(new_line[0], new_line[2], new_line[5], new_line[12], (new_line[15], new_line[16]), (new_line[17], new_line[18]))])
                        alin += 1
                        target_id = new_line[0]
                    else:
                        d[query_id] = [[(new_line[0], new_line[2], new_line[5], new_line[12], (new_line[15], new_line[16]), (new_line[17], new_line[18]))]]
                        alin = 0
                        target_id = new_line[0]
                else:
                    query_id = new_line[3]
                    line -= 1
            line += 1
        return d


    def readRSEMGene(self, expr_file):
        r = open(expr_file, "r")
        dic = {}
        r.readline()
        line = r.readline()
        while line:
            new_line = line.split("\t")
            if float(new_line[6]) > 0:
                dic[new_line[0]] = 0
            line = r.readline()
        r.close()
        return dic
    

    def readRSEMIso(self, expr_file):
        r = open(expr_file, "r")
        dic = {}
        r.readline()
        line = r.readline()
        while line:
            new_line = line.split("\t")
            dic[new_line[0]] = new_line[1].strip()
            line = r.readline()
        r.close()
        return dic
    

    def readCuffGene(self, expr_file):
        r = open(expr_file, "r")
        dic = {}
        r.readline()
        line = r.readline()
        while line:
            new_line = line.split("\t")
            if float(new_line[9]) > 0:
                dic[new_line[0]] = 0
            line = r.readline()
        r.close()
        return dic
    
    def readCuffIso(self, expr_file):
        r = open(expr_file, "r")
        dic = {}
        r.readline()
        line = r.readline()
        while line:
            new_line = line.split("\t")
            dic[new_line[0]] = new_line[3].strip()
            line = r.readline()
        r.close()
        return dic
    

    def run(self, ref_hmms):
        self.check_genes = []
        cegs_dict = {}
        for hmm in ref_hmms:
            if not self.hmmsearch_output_ref_as.has_key(hmm[0]):
                cegs_dict[hmm[0]] = ["-", "-", "-", "-", "-", "-"]
            else:
                result = self.identifyCeg(hmm[0])
                if result[4] == []:
                    cegs_dict[hmm[0]] = ["-", "-", "-", "-", "-", "-"]
                else:
                    self.check_genes.append(result[4][0])
                    cegs_dict[hmm[0]] = [hmm[1], self.hmmsearch_output_ref_as[hmm[0]][0][0][2], result[0], result[1], result[2], result[3]]
        return cegs_dict
    

    def identifyCeg(self, hmm):
        ceg_cov = 0.0
        protein_cov = 0.0
        length_protein = 0.0
        protein_id = ""
        gene_l = []
        for protein in self.hmmsearch_output_ref_as[hmm]:
            isoform = protein[0][0].split("|m")[0]
            gene = self.isoforms[isoform]
            if not self.genes.has_key(gene):
                continue
            if gene not in self.check_genes:
                domain = sorted(protein, key = lambda x: float(x[3]), reverse = False)[0]
                if float(int(domain[5][1]) - (int(domain[5][0]) - 1)) / int(domain[1]) >= self.prot_min_cov:
                    if float(int(domain[4][1]) - (int(domain[4][0]) - 1)) / int(domain[2]) >= self.ceg_min_cov:
                        gene_l.append(gene)
                        ceg_cov = float(int(domain[4][1]) - (int(domain[4][0]) - 1)) / int(domain[2])
                        protein_cov = float(int(domain[5][1]) - (int(domain[5][0]) - 1)) / int(domain[1])
                        length_protein = float(domain[1])
                        protein_id = domain[0]
                        return ceg_cov, protein_cov, length_protein, protein_id, gene_l
        return ceg_cov, protein_cov, length_protein, protein_id, gene_l


CEGs_con = [('KOG0002', '4'), ('KOG0018', '1'), ('KOG0025', '1'), ('KOG0062', '2'), ('KOG0077', '4'), ('KOG0122', '1'), ('KOG0142', '2'), ('KOG0174', '3'),
        ('KOG0175', '4'), ('KOG0176', '4'), ('KOG0177', '1'), ('KOG0179', '2'), ('KOG0180', '3'), ('KOG0181', '4'), ('KOG0182', '3'), ('KOG0184', '2'),
        ('KOG0188', '3'), ('KOG0209', '2'), ('KOG0211', '3'), ('KOG0225', '3'), ('KOG0233', '3'), ('KOG0261', '2'), ('KOG0271', '3'), ('KOG0276', '3'),
        ('KOG0279', '4'), ('KOG0285', '3'), ('KOG0291', '1'), ('KOG0292', '2'), ('KOG0302', '1'), ('KOG0329', '4'), ('KOG0346', '1'), ('KOG0357', '4'),
        ('KOG0358', '4'), ('KOG0359', '4'), ('KOG0361', '4'), ('KOG0363', '4'), ('KOG0364', '4'), ('KOG0365', '1'), ('KOG0366', '3'), ('KOG0376', '2'),
        ('KOG0400', '4'), ('KOG0419', '4'), ('KOG0424', '4'), ('KOG0434', '3'), ('KOG0462', '2'), ('KOG0466', '4'), ('KOG0469', '4'), ('KOG0477', '3'),
        ('KOG0481', '3'), ('KOG0524', '4'), ('KOG0556', '3'), ('KOG0559', '2'), ('KOG0563', '3'), ('KOG0567', '2'), ('KOG0650', '1'), ('KOG0687', '2'),
        ('KOG0688', '4'), ('KOG0729', '4'), ('KOG0741', '2'), ('KOG0780', '3'), ('KOG0815', '3'), ('KOG0820', '4'), ('KOG0861', '3'), ('KOG0862', '1'),
        ('KOG0871', '1'), ('KOG0876', '3'), ('KOG0888', '4'), ('KOG0894', '3'), ('KOG0927', '3'), ('KOG0937', '4'), ('KOG0948', '3'), ('KOG0960', '2'),
        ('KOG0964', '1'), ('KOG0969', '3'), ('KOG0985', '3'), ('KOG0989', '2'), ('KOG0991', '4'), ('KOG1058', '2'), ('KOG1068', '1'), ('KOG1088', '1'),
        ('KOG1099', '3'), ('KOG1112', '4'), ('KOG1123', '4'), ('KOG1137', '2'), ('KOG1145', '1'), ('KOG1159', '1'), ('KOG1185', '1'), ('KOG1211', '1'),
        ('KOG1235', '1'), ('KOG1241', '1'), ('KOG1272', '1'), ('KOG1299', '1'), ('KOG1322', '4'), ('KOG1335', '4'), ('KOG1349', '4'), ('KOG1350', '4'),
        ('KOG1353', '4'), ('KOG1355', '3'), ('KOG1358', '1'), ('KOG1367', '4'), ('KOG1373', '4'), ('KOG1374', '3'), ('KOG1390', '3'), ('KOG1393', '2'),
        ('KOG1394', '2'), ('KOG1439', '4'), ('KOG1458', '3'), ('KOG1463', '3'), ('KOG1466', '1'), ('KOG1468', '2'), ('KOG1498', '1'), ('KOG1523', '2'),
        ('KOG1532', '2'), ('KOG1533', '2'), ('KOG1534', '2'), ('KOG1535', '1'), ('KOG1540', '2'), ('KOG1541', '3'), ('KOG1549', '4'), ('KOG1555', '4'),
        ('KOG1556', '3'), ('KOG1596', '4'), ('KOG1597', '2'), ('KOG1636', '2'), ('KOG1637', '4'), ('KOG1643', '3'), ('KOG1647', '3'), ('KOG1662', '1'),
        ('KOG1664', '1'), ('KOG1688', '3'), ('KOG1712', '2'), ('KOG1723', '3'), ('KOG1727', '1'), ('KOG1733', '2'), ('KOG1746', '2'), ('KOG1758', '1'),
        ('KOG1760', '1'), ('KOG1769', '3'), ('KOG1770', '4'), ('KOG1774', '4'), ('KOG1775', '4'), ('KOG1780', '3'), ('KOG1781', '3'), ('KOG1795', '4'),
        ('KOG1800', '1'), ('KOG1816', '1'), ('KOG1872', '1'), ('KOG1885', '4'), ('KOG1889', '1'), ('KOG1936', '3'), ('KOG1942', '4'), ('KOG1980', '1'),
        ('KOG2004', '2'), ('KOG2017', '1'), ('KOG2035', '2'), ('KOG2036', '3'), ('KOG2044', '1'), ('KOG2104', '2'), ('KOG2303', '3'), ('KOG2311', '2'),
        ('KOG2415', '3'), ('KOG2446', '4'), ('KOG2451', '3'), ('KOG2472', '2'), ('KOG2481', '1'), ('KOG2519', '3'), ('KOG2529', '4'), ('KOG2531', '1'),
        ('KOG2535', '4'), ('KOG2537', '2'), ('KOG2572', '3'), ('KOG2575', '1'), ('KOG2606', '1'), ('KOG2613', '2'), ('KOG2623', '1'), ('KOG2638', '3'),
        ('KOG2653', '3'), ('KOG2680', '4'), ('KOG2703', '2'), ('KOG2707', '1'), ('KOG2719', '1'), ('KOG2726', '1'), ('KOG2728', '1'), ('KOG2757', '1'),
        ('KOG2770', '2'), ('KOG2775', '4'), ('KOG2781', '2'), ('KOG2784', '3'), ('KOG2785', '1'), ('KOG2792', '1'), ('KOG2807', '1'), ('KOG2825', '3'),
        ('KOG2833', '2'), ('KOG2851', '1'), ('KOG2874', '4'), ('KOG2906', '2'), ('KOG2908', '1'), ('KOG2909', '1'), ('KOG2916', '3'), ('KOG2930', '4'),
        ('KOG2948', '2'), ('KOG2967', '1'), ('KOG2971', '2'), ('KOG3013', '2'), ('KOG3048', '1'), ('KOG3049', '4'), ('KOG3052', '3'), ('KOG3147', '1'),
        ('KOG3157', '2'), ('KOG3163', '4'), ('KOG3167', '2'), ('KOG3174', '3'), ('KOG3180', '3'), ('KOG3185', '4'), ('KOG3189', '3'), ('KOG3205', '1'),
        ('KOG3218', '2'), ('KOG3222', '2'), ('KOG3232', '2'), ('KOG3237', '1'), ('KOG3239', '1'), ('KOG3275', '2'), ('KOG3285', '2'), ('KOG3295', '3'),
        ('KOG3297', '1'), ('KOG3313', '1'), ('KOG3318', '1'), ('KOG3330', '2'), ('KOG3343', '1'), ('KOG3349', '1'), ('KOG3361', '4'), ('KOG3387', '4'),
        ('KOG3400', '2'), ('KOG3404', '4'), ('KOG3405', '3'), ('KOG3432', '3'), ('KOG3448', '4'), ('KOG3459', '4'), ('KOG3479', '2'), ('KOG3482', '4'),
        ('KOG3493', '4'), ('KOG3497', '4'), ('KOG3499', '4'), ('KOG3503', '4'), ('KOG3855', '1'), ('KOG3954', '3'), ('KOG4392', '3'), ('KOG4655', '2')]


def writeTofile(result):
    f = open("cegs_identified.txt", "w")
    f.write("CEG\tConservation_level\tCEG_length\tCEG_coverage\tProtein_ID\tProtein_length\tProtein_coverage\n")
    keys = sorted(result.keys(), key = lambda x: int(x[3:]), reverse=False)
    identified = 0
    con1 = 0
    con2 = 0
    con3 = 0
    con4 = 0
    for key in keys:
        if str(result[key][0]) != "-":
            identified += 1
            if result[key][0] == "1":
                con1 += 1
            elif result[key][0] == "2":
                con2 += 1
            elif result[key][0] == "3":
                con3 += 1
            elif result[key][0] == "4":
                con4 += 1
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(key), result[key][0], result[key][1], str(result[key][2]), result[key][5], str(result[key][4]), str(result[key][3])))
    print "Identified: %s" % str(identified)
    print "Group 1: %s" % str(con1)
    print "Group 2: %s" % str(con2)
    print "Group 3: %s" % str(con3)
    print "Group 4: %s" % str(con4)
    #print str(identified)+"\t"+str(con1)+"\t"+str(con2)+"\t"+str(con3)+"\t"+str(con4)
    f.close()


def helpMessage():
    print "usage: seqQIidentifyCEGs.py --strategy <strategy_type> --hmmsearch_output <hmmsearch_output> --expression_iso <expression_file> --expression_gene <expression_file> --ceg_min_cov <#> --prot_min_cov <#>"
    print "-s --strategy\tStrategy type: 'reference-based' or 'denovo'"
    print "-H --hmmsearch_output\thmmsearch domtblout output: 248 CEGs -> putative proteins obtained by TransDecoder"
    print "-i --expression_iso\tFile containing the expression levels per isoform: from RSEM or Cufflinks"
    print "-g --expression_gene\tFile containing the expression levels per gene: from RSEM or Cufflinks"
    print "-p --ceg_min_cov\tMinimum coverage for profile HMMs: real number between 0-1 (0.8 by default)"
    print "-P --prot_min_cov\tMinimum coverage for protein sequences: real number between 0-1 (0.8 by default)"
    print "-h --help\tThis help message"
    sys.exit()

if __name__ == "__main__":
    
    start_time = time.time()
    
    inCommandLine = False
    
    if inCommandLine:
        if len(sys.argv) == 1: helpMessage()
        
        ceg_min_cov = 0.8
        prot_min_cov = 0.8
        
        try:
            opts, args = getopt.getopt(sys.argv[1:],"hs:H:i:g:p:P:", ["help", "strategy=", "hmmsearch_output_ref_as=", "expression_iso=", "expression_gene=", "ceg_min_cov=", "prot_min_cov="])
        except getopt.GetoptError:
            helpMessage()
            
        for opt, arg in opts:
            if opt in ('-h', "--help"): helpMessage()
            elif opt in ("-s", "--strategy"): strategy = arg
            elif opt in ("-H", "--hmmsearch_output_ref_as"): hmmsearch_output_ref_as = arg
            elif opt in ("-i", "--expression_iso"): expression_iso = arg
            elif opt in ("-g", "--expression_gene"): expression_gene = arg
            elif opt in ("-p", "--ceg_min_cov"): ceg_min_cov = float(arg)
            elif opt in ("-P", "--prot_min_cov"): prot_min_cov = float(arg)
    else:
        strategy = "reference-based"
        hmmsearch_output_ref_as = "/Volumes/KINGSTON/Transdecoder_SS_hmmer_ref_100_SC/hmmsearch_ref_hmms_assem_trpts.txt"
        expression_iso = "/Volumes/KINGSTON/Transdecoder_SS_hmmer_ref_100_SC/isoforms.fpkm_tracking"
        expression_gene = "/Volumes/KINGSTON/Transdecoder_SS_hmmer_ref_100_SC/genes.fpkm_tracking"
        ceg_min_cov = 0.8
        prot_min_cov = 0.8
    
    try:
        cegs = MyCegs(strategy, hmmsearch_output_ref_as, expression_iso, expression_gene, ceg_min_cov, prot_min_cov)
    except (NameError, IOError):
        helpMessage()
    
    result = cegs.run(CEGs_con)
    writeTofile(result)
    
    print "\nFinish!"
    
    print "\nA tab delimited text file was written to cegs_identified.txt"
