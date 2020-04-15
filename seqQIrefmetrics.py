'''
@author: abelsousa
'''

import sys, getopt, time
from Bio.Blast import NCBIXML

class MyMetrics:
    
    def __init__(self, strategy, blast_output_rf_as, blast_output_as_rf, expression_file, assem_min_align, cont_ref_min_align, chim_min_align):
        self.blast_output_rf_as = blast_output_rf_as
        self.blast_output_as_rf = blast_output_as_rf
        try:
            self.assem_min_align = float(assem_min_align)
        except ValueError:
            print "Unrecognized value! Please input a real number between 0-1"
            sys.exit()
        try:
            self.cont_ref_min_align = float(cont_ref_min_align)
        except ValueError:
            print "Unrecognized value! Please input a real number between 0-1"
            sys.exit()
        try:
            self.chim_min_align = chim_min_align
        except ValueError:
            print "Unrecognized value! Please input a real number between 0-1"
            sys.exit()
        self.ref, self.assem = self.calculateRefAssemTrpts()
        if strategy == "reference-based":
            self.isoforms = self.readCuffIso(expression_file)
            self.express = len(self.isoforms)
        elif strategy == "denovo":
            self.isoforms = self.readRSEMIso(expression_file)
            self.express = len(self.isoforms)
        else:
            print "This strategy type is not allowed! Please input 'reference-based' or 'denovo' (without quotes)"
            sys.exit()


    def calculateRefAssemTrpts(self):
        try:
            handle_rf_as = open(self.blast_output_rf_as)
        except TypeError:
            print "Unrecognized file! Please input a xml file"
            sys.exit()
        try:
            handle_as_rf = open(self.blast_output_as_rf)
        except TypeError:
            print "Unrecognized file! Please input a xml file"
            sys.exit()
        try:
            n_rf = next(NCBIXML.parse(handle_as_rf)).database_sequences
        except ValueError:
            print "Unrecognized file! Please input a xml file"
            sys.exit()
        try:
            n_as = next(NCBIXML.parse(handle_rf_as)).database_sequences
        except ValueError:
            print "Unrecognized file! Please input a xml file"
            sys.exit()
        handle_rf_as.close(), handle_as_rf.close()
        return n_rf, n_as
    

    def readCuffIso(self, expr_file):
        try:
            r = open(expr_file, "r")
        except TypeError:
            print "Unrecognized file! Please input a Cufflinks expression file per isoform"
            sys.exit()
        expressedIso = {}
        header = r.readline()
        if len(header.split("\t")) != 13:
            print "Unrecognized file! Please input a Cufflinks expression file per isoform"
            sys.exit()
        line = r.readline()
        while line:
            if float(line.split("\t")[9]) > 0:
                expressedIso[line.split("\t")[0]] = float(line.split("\t")[7]), float(line.split("\t")[9]), "-"
            line = r.readline()
        r.close()
        return expressedIso
    

    def readRSEMIso(self, expr_file):
        try:
            r = open(expr_file, "r")
        except TypeError:
            print "Unrecognized file! Please input a RSEM expression file per isoform"
            sys.exit()
        expressedIso = {}
        header = r.readline()
        if len(header.split("\t")) != 8:
            print "Unrecognized file! Please input a RSEM expression file per isoform"
            sys.exit()
        line = r.readline()
        while line:
            if float(line.split("\t")[6]) > 0:
                expressedIso[line.split("\t")[0]] = float(line.split("\t")[2]), float(line.split("\t")[6]), float(line.split("\t")[4])
            line = r.readline()
        r.close()
        return expressedIso
    
    
    def runMetrics(self):
        cont = open("contiguity.txt", "w")
        cont.write("ID\tLength\tFPKM\tExpected counts\n")
        frag1 = open("fragmentation_1.txt", "w")
        frag1.write("ID\tLength\tFPKM\tExpected counts\n")
        frag2 = open("fragmentation_2.txt", "w")
        frag2.write("ID\tLength\tFPKM\tExpected counts\n")
        frag3 = open("fragmentation_3.txt", "w")
        frag3.write("ID\tLength\tFPKM\tExpected counts\n")
        frag4 = open("fragmentation_4.txt", "w")
        frag4.write("ID\tLength\tFPKM\tExpected counts\n")
        frag5 = open("fragmentation_5.txt", "w")
        frag5.write("ID\tLength\tFPKM\tExpected counts\n")
        handle_rf_as = open(self.blast_output_rf_as)
        blast_records = NCBIXML.parse(handle_rf_as)
        cov_sum = 0.0
        align_length_sum = 0.0
        corr_bases_sum = 0.0
        hits_ref_sum = 0.0
        cont_sum = 0.0
        frag_sum = 0.0
        frag_sum_1 = 0.0
        frag_sum_2 = 0.0
        frag_sum_3 = 0.0
        frag_sum_4 = 0.0
        frag_sum_5 = 0.0
        assem_trpts_used = 0.0
        assem_trpts_used_frag = 0.0
        self.check_trpts = {}
        for blast_record in blast_records:
            result_iter = self.identifiedMetrics(blast_record)
            cov_sum += result_iter[0]
            align_length_sum += result_iter[1]
            corr_bases_sum += result_iter[2]
            if result_iter[3] == 1:
                hits_ref_sum += 1
                cont_sum += 1
                assem_trpts_used += 1
                cont.write("%s\t%s\t%s\t%s\n" % (result_iter[5][0], str(self.isoforms[result_iter[5][0]][0]), str(self.isoforms[result_iter[5][0]][1]), str(self.isoforms[result_iter[5][0]][2])))
            elif result_iter[4] == 1:
                hits_ref_sum += 1
                frag_sum += 1
                frag_sum_1 += 1
                assem_trpts_used += result_iter[4]
                assem_trpts_used_frag += result_iter[4]
                for i in result_iter[5]:
                    frag1.write("%s\t%s\t%s\t%s\n" % (i, str(self.isoforms[i][0]), str(self.isoforms[i][1]), str(self.isoforms[i][2])))
            elif result_iter[4] == 2:
                hits_ref_sum += 1
                frag_sum += 1
                frag_sum_2 += 1
                assem_trpts_used += result_iter[4]
                assem_trpts_used_frag += result_iter[4]
                for i in result_iter[5]:
                    frag2.write("%s\t%s\t%s\t%s\n" % (i, str(self.isoforms[i][0]), str(self.isoforms[i][1]), str(self.isoforms[i][2])))
            elif result_iter[4] == 3:
                hits_ref_sum += 1
                frag_sum += 1
                frag_sum_3 += 1
                assem_trpts_used += result_iter[4]
                assem_trpts_used_frag += result_iter[4]
                for i in result_iter[5]:
                    frag3.write("%s\t%s\t%s\t%s\n" % (i, str(self.isoforms[i][0]), str(self.isoforms[i][1]), str(self.isoforms[i][2])))
            elif result_iter[4] == 4:
                hits_ref_sum += 1
                frag_sum += 1
                frag_sum_4 += 1
                assem_trpts_used += result_iter[4]
                assem_trpts_used_frag += result_iter[4]
                for i in result_iter[5]:
                    frag4.write("%s\t%s\t%s\t%s\n" % (i, str(self.isoforms[i][0]), str(self.isoforms[i][1]), str(self.isoforms[i][2])))
            elif result_iter[4] >= 5:
                hits_ref_sum += 1
                frag_sum += 1
                frag_sum_5 += 1
                assem_trpts_used += result_iter[4]
                assem_trpts_used_frag += result_iter[4]
                for i in result_iter[5]:
                    frag5.write("%s\t%s\t%s\t%s\n" % (i, str(self.isoforms[i][0]), str(self.isoforms[i][1]), str(self.isoforms[i][2])))
        handle_rf_as.close(), cont.close(), frag1.close(), frag2.close(), frag3.close(), frag4.close(), frag5.close()
        identified = 100 * (hits_ref_sum / self.ref) 
        completeness = 100 * (cov_sum / hits_ref_sum)
        contiguity = 100 * (cont_sum / hits_ref_sum)
        fragmented = 100 * (frag_sum / hits_ref_sum)
        fragmented_1 = 100 * (frag_sum_1 / hits_ref_sum)
        fragmented_2 = 100 * (frag_sum_2 / hits_ref_sum)
        fragmented_3 = 100 * (frag_sum_3 / hits_ref_sum)
        fragmented_4 = 100 * (frag_sum_4 / hits_ref_sum)
        fragmented_5 = 100 * (frag_sum_5 / hits_ref_sum)
        accuracy = 100 * (corr_bases_sum / align_length_sum)  
        result_completenessCont = identified, hits_ref_sum, completeness, cov_sum, contiguity, cont_sum, fragmented, frag_sum, assem_trpts_used, assem_trpts_used_frag, accuracy, [fragmented_1, frag_sum_1, fragmented_2, frag_sum_2, fragmented_3, frag_sum_3, fragmented_4, frag_sum_4, fragmented_5, frag_sum_5]
        handle_as_rf = open(self.blast_output_as_rf)
        non_mat = open("non_match.txt", "w")
        non_mat.write("ID\tLength\tFPKM\tExpected counts\n")
        chim = open("chimerism.txt", "w")
        chim.write("ID\tLength\tFPKM\tExpected counts\n")
        blast_records = NCBIXML.parse(handle_as_rf)
        chimaeras = 0.0
        no_hits = 0.0
        
        for blast_record in blast_records:
            result_iter = self.chimerismNonMatch(blast_record)
            if result_iter[0] == 1:
                chimaeras += 1
                chim.write("%s\t%s\t%s\t%s\n" % (str(blast_record.query).split(" ")[0], str(self.isoforms[str(blast_record.query).split(" ")[0]][0]), str(self.isoforms[str(blast_record.query).split(" ")[0]][1]), str(self.isoforms[str(blast_record.query).split(" ")[0]][2])))
            elif result_iter[1] == 1:
                no_hits += 1
                non_mat.write("%s\t%s\t%s\t%s\n" % (str(blast_record.query).split(" ")[0], str(self.isoforms[str(blast_record.query).split(" ")[0]][0]), str(self.isoforms[str(blast_record.query).split(" ")[0]][1]), str(self.isoforms[str(blast_record.query).split(" ")[0]][2])))
        handle_as_rf.close(), non_mat.close(), chim.close()
        perc_chimaeras = (chimaeras / float(self.express)) * 100
        perc_no_hits = (no_hits / float(self.express)) * 100
        result_chimerism_ord = perc_chimaeras, chimaeras, perc_no_hits, no_hits
        return result_completenessCont, result_chimerism_ord
    

    def identifiedMetrics(self, record):
        align_length = 0.0
        corr_bases = 0.0
        frag_hits = 0.0
        cov = 0.0
        cont = 0.0
        aligns = []
        trpts = []
        for alignment in xrange(len(record.alignments)):
            if not self.isoforms.has_key(str(record.alignments[alignment].hit_def).split(" ")[0]):
                continue
            if not self.check_trpts.has_key(str(record.alignments[alignment].hit_def)):
                if record.alignments[alignment].hsps[0].align_length >= record.alignments[alignment].length:
                    p_align_assem = 1.0
                else:
                    p_align_assem = float(record.alignments[alignment].hsps[0].align_length) / record.alignments[alignment].length
                if p_align_assem >= self.assem_min_align:
                    if record.alignments[alignment].hsps[0].align_length >= record.query_length:
                        p_align_ref = 1.0
                    else:
                        p_align_ref = float(record.alignments[alignment].hsps[0].align_length) / record.query_length
                    if p_align_ref >= self.cont_ref_min_align:
                        cov += p_align_ref
                        cont += 1.0
                        align_length += record.alignments[alignment].hsps[0].align_length
                        corr_bases += record.alignments[alignment].hsps[0].identities
                        self.check_trpts[str(record.alignments[alignment].hit_def)] = 0
                        trpts.append(str(record.alignments[alignment].hit_def).split(" ")[0])
                        return cov, align_length, corr_bases, cont, frag_hits, trpts
                    else:
                        aligns.append((alignment, record.alignments[alignment].hsps[0].query_start))
        if len(aligns) == 0:
            return cov, align_length, corr_bases, cont, frag_hits, trpts
        aligns = sorted(aligns, key = lambda x: x[1], reverse = False)
        self.check_trpts[str(record.alignments[aligns[0][0]].hit_def)] = 0
        trpts.append(str(record.alignments[aligns[0][0]].hit_def).split(" ")[0])
        align_length += record.alignments[aligns[0][0]].hsps[0].align_length
        corr_bases += record.alignments[aligns[0][0]].hsps[0].identities
        frag_hits += 1.0
        cov += float(record.alignments[aligns[0][0]].hsps[0].align_length) / record.query_length
        end_region = record.alignments[aligns[0][0]].hsps[0].query_end
        for align in xrange(len(aligns)-1):
            if aligns[align+1][1] > end_region:
                self.check_trpts[str(record.alignments[aligns[align+1][0]].hit_def)] = 0
                trpts.append(str(record.alignments[aligns[align+1][0]].hit_def).split(" ")[0])
                align_length += record.alignments[aligns[align+1][0]].hsps[0].align_length
                corr_bases += record.alignments[aligns[align+1][0]].hsps[0].identities
                frag_hits += 1
                cov += float(record.alignments[aligns[align+1][0]].hsps[0].align_length) / record.query_length
                if cov >= 1.0:
                    cov = 1.0
                    return cov, align_length, corr_bases, cont, frag_hits, trpts
                end_region = record.alignments[aligns[align+1][0]].hsps[0].query_end
            elif not record.alignments[aligns[align+1][0]].hsps[0].query_end <= end_region:
                self.check_trpts[str(record.alignments[aligns[align+1][0]].hit_def)] = 0
                trpts.append(str(record.alignments[aligns[align+1][0]].hit_def).split(" ")[0])
                align_length += (record.alignments[aligns[align+1][0]].hsps[0].query_end - end_region)
                corr_bases += record.alignments[aligns[align+1][0]].hsps[0].match[record.alignments[aligns[align+1][0]].hsps[0].align_length - (record.alignments[aligns[align+1][0]].hsps[0].query_end - end_region):].count("|")
                frag_hits += 1
                cov += float(record.alignments[aligns[align+1][0]].hsps[0].query_end - end_region) / record.query_length
                if cov >= 1.0:
                    cov = 1.0
                    return cov, align_length, corr_bases, cont, frag_hits, trpts
                end_region = record.alignments[aligns[align+1][0]].hsps[0].query_end
            else:
                self.check_trpts[str(record.alignments[aligns[align+1][0]].hit_def)] = 0
                trpts.append(str(record.alignments[aligns[align+1][0]].hit_def).split(" ")[0])
                align_length += record.alignments[aligns[align+1][0]].hsps[0].align_length
                corr_bases += record.alignments[aligns[align+1][0]].hsps[0].identities
                frag_hits += 1
        return cov, align_length, corr_bases, cont, frag_hits, trpts


    def chimerismNonMatch(self, record):
        chimaera = 0
        no_hits = 0
        if not self.isoforms.has_key(str(record.query).split(" ")[0]):
            return chimaera, no_hits
        if self.check_trpts.has_key(str(record.query)):
            return chimaera, no_hits
        aligns = []
        if len(record.alignments) == 0:
            no_hits += 1
            return chimaera, no_hits
        for alignment in xrange(len(record.alignments)):
            if record.alignments[alignment].hsps[0].align_length >= record.alignments[alignment].length:
                p_align_ref = 1.0
            else:
                p_align_ref = float(record.alignments[alignment].hsps[0].align_length) / record.alignments[alignment].length
            if p_align_ref >= self.chim_min_align:
                aligns.append((alignment, record.alignments[alignment].hsps[0].query_start))
        if len(aligns) > 1:
            aligns = sorted(aligns, key = lambda x: x[1], reverse=False)
            end_region = record.alignments[aligns[0][0]].hsps[0].query_end
            for align in xrange(len(aligns)-1):
                if aligns[align+1][1] > end_region:
                    chimaera += 1
                    return chimaera, no_hits
                elif not record.alignments[aligns[align+1][0]].hsps[0].query_end <= end_region:
                    end_region = record.alignments[aligns[align+1][0]].hsps[0].query_end
        no_hits += 1
        return chimaera, no_hits


def helpMessage():
    print "usage: seqQIrefmetrics.py --strategy <strategy_type> --blast_output_rf_as <blast_output> --blast_output_as_rf <blast_output> --expression_file <expression_file> --assembl_min_cov <#> --ref_cont_min_cov <#> --ref_chim_min_cov <#>"
    print "-s --strategy\tStrategy type: 'reference-based' or 'denovo'"
    print "-r --blast_output_rf_as\tXML blast output: reference transcripts -> assembled transcripts"
    print "-a --blast_output_as_rf\tXML blast output: assembled transcripts -> reference transcripts"
    print "-e --expression_file\tFile containing the transcripts expression levels per isoform: from RSEM or Cufflinks"
    print "-A --assembl_min_cov\tMinimum coverage for assembled transcripts: real number between 0-1 (0.8 by default)"
    print "-C --ref_cont_min_cov\tMinimum coverage for reference transcripts for contiguity: real number between 0-1 (0.8 by default)"
    print "-c --ref_chim_min_cov\tMinimum coverage for reference transcripts for chimerism: real number between 0-1 (0.8 by default)"
    print "-h --help\tThis help message"
    sys.exit()

if __name__ == "__main__":
    
    start_time = time.time()
    
    inCommandLine = False
    
    if inCommandLine:
        if len(sys.argv) == 1: helpMessage()
        
        assembl_min_cov = 0.8
        ref_cont_min_cov = 0.8
        ref_chim_min_cov = 0.8
                
        try:
            opts, args = getopt.getopt(sys.argv[1:],"hs:r:a:e:A:C:c:", ["help", "strategy=",  "blast_output_rf_as=", "blast_output_as_rf=", "expression_file=", "assembl_min_cov=", "ref_cont_min_cov=", "ref_chim_min_cov="])
        except getopt.GetoptError:
            helpMessage()
            
        for opt, arg in opts:
            if opt in ('-h', "--help"): helpMessage()
            elif opt in ("-s", "--strategy"): strategy = arg
            elif opt in ("-r", "--blast_output_rf_as"): blast_output_rf_as = arg
            elif opt in ("-a", "--blast_output_as_rf"): blast_output_as_rf = arg
            elif opt in ("-e", "--expression_file"): expression_file = arg
            elif opt in ("-A", "--assembl_min_cov"): assembl_min_cov = arg
            elif opt in ("-C", "--ref_cont_min_cov"): ref_cont_min_cov = arg
            elif opt in ("-c", "--ref_chim_min_cov"): ref_chim_min_cov = arg
        
    else:
        strategy = "reference-based"
        blast_output_rf_as = "/Volumes/KINGSTON/BLAST_evalue_1e-6_ref_100_SC/blast_results_rf_as_refassembly_1.xml"
        blast_output_as_rf = "/Volumes/KINGSTON/BLAST_evalue_1e-6_ref_100_SC/blast_results_as_rf_refassembly_1.xml"
        expression_file = "/Volumes/KINGSTON/BLAST_evalue_1e-6_ref_100_SC/isoforms.fpkm_tracking"
        assembl_min_cov = 0.8
        ref_cont_min_cov = 0.8
        ref_chim_min_cov = 0.8
    
    try:
        metrics = MyMetrics(strategy, blast_output_rf_as, blast_output_as_rf, expression_file, assembl_min_cov, ref_cont_min_cov, ref_chim_min_cov)
    except (NameError, IOError):
        helpMessage()
    print "Number of reference transcripts: %s" % (str(metrics.ref))
    print "Number of assembled transcripts: %s" % (str(metrics.assem))
    print "Number of expressed transcripts: %s" % (str(metrics.express))
    print
    print "Minimum coverage for assembled transcripts: %s" % (str(metrics.assem_min_align))
    print "Minimum coverage for reference transcripts (contiguity): %s" % (str(metrics.cont_ref_min_align))
    print "Minimum coverage for reference transcripts (chimerism): %s" % (str(metrics.chim_min_align))
    print
    print "Calculating..."
    
    result = metrics.runMetrics()
    
    print
    print "Identification -> %s%% (%s reference transcripts)" % (str(result[0][0]), str(result[0][1]))
    print "Coverage -> %s%%" % (str(result[0][2]))
    print "Contiguity -> %s%% (%s reference transcripts)" % (str(result[0][4]), str(result[0][5]))
    print "Fragmentation -> %s%% (%s reference transcripts)" % (str(result[0][6]), str(result[0][7]))
    print "Fragmentation(1) -> %s%% (%s reference transcripts)" % (str(result[0][11][0]), str(result[0][11][1]))
    print "Fragmentation(2) -> %s%% (%s reference transcripts)" % (str(result[0][11][2]), str(result[0][11][3]))
    print "Fragmentation(3) -> %s%% (%s reference transcripts)" % (str(result[0][11][4]), str(result[0][11][5]))
    print "Fragmentation(4) -> %s%% (%s reference transcripts)" % (str(result[0][11][6]), str(result[0][11][7]))
    print "Fragmentation(5) -> %s%% (%s reference transcripts)" % (str(result[0][11][8]), str(result[0][11][9]))
    print "Accuracy -> %s%%" % (str(result[0][10]))
    print "Chimerism -> %s%% (%s assembled transcripts)" % (str(result[1][0]), str(result[1][1]))
    print "Non-match -> %s%% (%s assembled transcripts)" % (str(result[1][2]), str(result[1][3]))
    print
    print "Total number of assembled transcripts used for fragmentation -> %s" % (str(result[0][9]))
    print "Total number of assembled transcripts used for contiguity/fragmentation -> %s" % (str(result[0][8]))
    print
    #print str(metrics.assem)+"\t"+str(metrics.express)+"\t"+str(result[0][1])+"\t"+str(result[0][0])+"\t"+str(result[0][2])+"\t"+str(result[0][5])+"\t"+str(result[0][4])+"\t"+str(result[0][7])+"\t"+str(result[0][6])+"\t"+str(result[0][10])+"\t"+str(result[1][1])+"\t"+str(result[1][0])+"\t"+str(result[1][3])+"\t"+str(result[1][2])+"\t"+str(result[0][9])+"\t"+str(result[0][8])
    
    # write results to a file
    f = open("reports.txt", "w")
    f.write("Number of reference transcripts: %s\n" % (str(metrics.ref)))
    f.write("Number of assembled transcripts: %s\n" % (str(metrics.assem)))
    f.write("Number of expressed transcripts: %s\n" % (str(metrics.express)))
    
    f.write("\nMinimum coverage for assembled transcripts: %s\n" % (str(metrics.assem_min_align)))
    f.write("Minimum coverage for reference transcripts (contiguity): %s\n" % (str(metrics.cont_ref_min_align)))
    f.write("Minimum coverage for reference transcripts (chimerism): %s\n" % (str(metrics.chim_min_align)))
    
    f.write("Identification -> %s%% (%s transcripts)\n" % (str(result[0][0]), str(result[0][1])))
    f.write("Coverage -> %s%%\n" % (str(result[0][2])))
    f.write("Contiguity -> %s%% (%s transcripts)\n" % (str(result[0][4]), str(result[0][5])))
    f.write("Fragmentation -> %s%% (%s transcripts)\n" % (str(result[0][6]), str(result[0][7])))
    f.write("Fragmentation(1) -> %s%% (%s transcripts)\n" % (str(result[0][11][0]), str(result[0][11][1])))
    f.write("Fragmentation(2) -> %s%% (%s transcripts)\n" % (str(result[0][11][2]), str(result[0][11][3])))
    f.write("Fragmentation(3) -> %s%% (%s transcripts)\n" % (str(result[0][11][4]), str(result[0][11][5])))
    f.write("Fragmentation(4) -> %s%% (%s transcripts)\n" % (str(result[0][11][6]), str(result[0][11][7])))
    f.write("Fragmentation(5) -> %s%% (%s transcripts)\n" % (str(result[0][11][8]), str(result[0][11][9])))
    f.write("Accuracy -> %s%%\n" % (str(result[0][10])))
    f.write("Chimerism -> %s%% (%s assembled transcripts)\n" % (str(result[1][0]), str(result[1][1])))
    f.write("Non-match -> %s%% (%s assembled transcripts)\n" % (str(result[1][2]), str(result[1][3])))
    f.write("\nTotal number of assembled transcripts used for fragmentation -> %s\n" % (str(result[0][9])))
    f.write("Total number of assembled transcripts used during contiguity/fragmentation -> %s\n" % (str(result[0][8])))
    f.write("\n"+str(metrics.assem)+"\t"+str(metrics.express)+"\t"+str(result[0][1])+"\t"+str(result[0][0])+"\t"+str(result[0][2])+"\t"+str(result[0][5])+"\t"+str(result[0][4])+"\t"+str(result[0][7])+"\t"+str(result[0][6])+"\t"+str(result[0][10])+"\t"+str(result[1][1])+"\t"+str(result[1][0])+"\t"+str(result[1][3])+"\t"+str(result[1][2])+"\t"+str(result[0][9])+"\t"+str(result[0][8]))
    f.close()

    seconds = time.time() - start_time
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    
    print "\nFinish!"
    print "%d:%02d:%02d" % (h, m, s)
    
    print "\nMetrics results saved to reports.txt"
    print "\nTranscripts length, FPKM and counts also saved to 'metric'.txt"