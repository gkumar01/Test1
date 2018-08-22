import os
import sys
import warnings
import re
import pysam
import pybedtools as bd
import numpy as np
import pandas as pd
from optparse import OptionParser

#https://gist.github.com/amblina/5ed9a61ce74ad90668f4e29d62e7eb79
class BamFilter(object):
    """ BAM File filter class """
   
    def __init__(self,
                 input_bam_fn=None,
                 target_region_fn = None
                 ):

        self.input_bam_fn = input_bam_fn
        self.target_region_fn = target_region_fn

        self.__minimum_base_quality = 30
        self.__minimum_mapq = 20
        self.__minimum_read_quality = 20

        print ("xx",self.input_bam_fn)
        print ("xx",self.target_region_fn)


    def get_avg_seq_qual_n_aln_read_count(self,bamfile):

        """ parse bam file to print % average sequencing quality
         Args: pysam bamfile
         Returns: avg_seq_qual = float, aligned_read_cnt = interger
        """

        qual_values = []
        aligned_read_cnt = 0

        for rd in bamfile.fetch(until_eof=True):
            rd_name = rd.qname
            rd_qual = np.mean(rd.query_alignment_qualities)
            #print(rd_name,rd_qual)
            qual_values.append(rd_qual)

            #check if read is mapped to reference genome
            if not rd.is_unmapped:
                aligned_read_cnt += 1


        avg_seq_qual = np.mean(qual_values)

        # print("1.Average sequencing Quality:%2.2f%%" %(avg_seq_qual))
        # print(" Total Aligned Read:%d" %(aligned_read_cnt))
        del qual_values[:]

        return (avg_seq_qual,aligned_read_cnt)


    # def aligned_read_cnt(self,bamfile):
    #
    #     """ count and display the number of aligned reads to reference sequence """
    #     aligned_read_cnt = 0
    #
    #     for rd in bamfile.fetch(until_eof=True):
    #         # rd_name = rd.qname
    #         # print("query:%s" % (rd_name))
    #
    #         if not rd.is_unmapped:
    #             aligned_read_cnt += 1
    #
    #     print(" Total Aligned Read:%d" %(aligned_read_cnt))
    #     return aligned_read_cnt

    def get_read_mapped_loc(self, rd, header):
        """ calculate read mapped coordinates """
        start = rd.pos + 1
        end = rd.pos + len(rd.seq)
        chrom = header['SQ'][rd.tid]['SN']
        location = (chrom, start, end)
        return location

    def get_pct_read_enrich_target(self,bamfile,mapped_read_cnt=0.000001):
        """ function to get enrichment of target sites
        Args: bamfile = read aligned bamfile,
            mapped_read_cnt = integer value containing number of aligned read

        Returns: float, percentage of aligned reads on the target
        """
        target = bd.BedTool(self.target_region_fn)
        target = target.remove_invalid().saveas()

        bam = bd.BedTool(self.input_bam_fn)
        bam = bam.bam_to_bed(stream=True)

        target = target.remove_invalid().saveas()
        tmp = bam.remove_invalid().saveas()

        ontarget_reads = tmp.coverage(target, counts=True).sort()
        target_read = np.array([int(ln[-1]) for ln in ontarget_reads])
        target_read_count = target_read.sum()

        #print(target_read_count)
        #mapped_read_cnt = self.aligned_read_cnt(bamfile)
        #print(mapped_read_cnt)
        #print(type(mapped_read_cnt))

        pct_ontarget = 100*(target_read_count/float(mapped_read_cnt))

        return pct_ontarget

    def get_uniform_coverage_target_base_cnt(self,bamfile,tbl):

        target_base_uniform_cov = 0
        for r_idx, row in tbl.iterrows():
            if r_idx == len(tbl):
                break
            chrom = tbl.iloc[r_idx]['Chromosome']
            start = int(tbl.iloc[r_idx]['Start'])
            end = int(tbl.iloc[r_idx]['Stop'])
            mean_cov = 0.2*float(tbl.iloc[r_idx]['MeanCoverage'])
            #print(mean_cov)
            for pileupcol in bamfile.pileup(chrom, start, end):
                pos, rd_num = pileupcol.pos, pileupcol.n
                #print ("\ncoverage at base %s = %s" %(pos,rd_num))
                if(rd_num < mean_cov ):
                    target_base_uniform_cov +=1

        #print(target_base_uniform_cov)
        return target_base_uniform_cov

    def get_target_summary_tbl(self,bamfile):
        target_summary_stat = []
        all_target_base_cnt = 0
        align_target_base_cnt = 0
        with open(self.target_region_fn,'r') as target:
            read_data = target.readlines()
            for ln in read_data:

                # print ln.strip('\n')
                location = ln.strip('\n').split('\t')
                chr,start,end = location[:3]

                #chr = re.sub('chr','',chr)
                start = int(start)
                end = int(end)

                #print(chr,start,end)
                cov_list = []
                read_cov_info=location

                align_base = 0
                try:
                    align_read = bamfile.fetch(chr,start,end).next()
                    align_base = len(align_read.positions)
                except StopIteration:
                    pass

                align_target_base_cnt += align_base
                target_len = end - start
                all_target_base_cnt += target_len

                #print(all_target_base_cnt, all_target_base_cnt)
                for pileupcol in bamfile.pileup(chr,start,end):
                    pos,rd_num = pileupcol.pos,pileupcol.n
                    cov_list.append(rd_num)

                #print(cov_list)

                mean_cov = np.mean(cov_list)
                mean_cov = "%.2f" % mean_cov
                sd_cov = np.std(cov_list)
                sd_cov = "%.2f" % sd_cov
                read_cov_info.extend([mean_cov,sd_cov])
                #print read_cov_info
                target_summary_stat.append(read_cov_info)

        percnt_align_target_base = 100*align_target_base_cnt/float(all_target_base_cnt)
        print("3.Percentage of base enrich in target:%2.2f%%" %(percnt_align_target_base))

        tbl = pd.DataFrame(target_summary_stat)
        tbl.columns=['Chromosome','Start','Stop','RegionID','MeanCoverage','StdDevCoverage']

        #print(tbl.head())


        uniform_base_cov_cnt = self.get_uniform_coverage_target_base_cnt(bamfile,tbl)
        uniform_base_cov_percnt = 100*uniform_base_cov_cnt/float(all_target_base_cnt)
        print("4.Percentage of Uniform coverage of base in target:%2.2f%%" % (uniform_base_cov_percnt))


        #mean_cov_allbase_target = tbl['MeanCoverage'].mean()
        #print("5. Mean Coverage for each base in all target region:",mean_cov_allbase_target)
        #print(tbl.head(5))
        #6. Create table
        tbl.to_csv("TEST1_target_summary_stat.tsv",sep="\t")

    def filter(self):

        if not os.path.exists(self.input_bam_fn) or \
                not os.path.exists(self.target_region_fn):
            print(" Check input files:")
            return

        bamfile = pysam.AlignmentFile(self.input_bam_fn, 'rb')

        (avg_seq_qual,aln_read_cnt) = self.get_avg_seq_qual_n_aln_read_count(bamfile)
        print("1.Average sequencing Quality:%2.2f%%" % (avg_seq_qual))

        pct_ontarget = self.get_pct_read_enrich_target(bamfile,aln_read_cnt)
        print("2. Percentage of read enriched in the target region:%2.2f%%" %(pct_ontarget))
        self.get_target_summary_tbl(bamfile)

        bamfile.close()

    def run(self):
        self.filter()


def get_params():
    parser = OptionParser()
    parser.add_option("-i", "--input-bam-fn",
                      dest="input_bam_fn",
                      default=None,
                      help="The input bam filename, eg. TEST1.bam")
    parser.add_option("-t", "--target-region-fn",
                      dest="target_region_fn",
                      default=None,
                      help="The input target file, TEST1_region.bed")

    (options, args) = parser.parse_args()
    return options, args

def run(options,args):
    bmf = BamFilter(input_bam_fn=options.input_bam_fn,
                    target_region_fn=options.target_region_fn)
    bmf.run()


if __name__ == '__main__':

    """
    python test_bam.py -i TEST1.bam -t TEST1_region.bed
    """
    options, args = get_params()
    run(options,args)
