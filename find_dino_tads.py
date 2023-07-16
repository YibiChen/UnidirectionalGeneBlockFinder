import sys
import os
from Bio import SeqIO
import matplotlib.pyplot as plt

class UniDirectionGeneBlock:
    """A unidirectional gene block is characterized by a cluster of genes that are all oriented in the same coding direction."""
    def __init__(self):
        self.genes = []
        self.start = None
        self.end = None
        self.scaffold = None
        self.orientation = None
    def __len__(self):
        return len(self.genes)
    def add_gene(self,gene_id,scaffold,start,end,orientation):
        if self.scaffold == None:
            self.scaffold = scaffold
        elif self.scaffold != scaffold:
            raise ValueError
        
        if self.orientation == None:
            self.orientation = orientation
        elif self.orientation != orientation:
            raise ValueError

        self.genes.append(gene_id)
        if self.start != None:
            self.start = min(self.start,start)
        else:
            self.start=start
        if self.end != None:
            self.end = max(self.end,end)
        else:
            self.end = end


def search_unidirectional_gene_blocks(gff3_file):
    """Idenfity all unidirectional gene blocks in the genome."""
    blocks = {}
    with open(gff3_file) as f:
        last_orientation = ""
        last_scaffold = ""
        new_block = UniDirectionGeneBlock()
        for line in f:
            record = line.strip().split()
            scaffold_id=record[0]
            orientation = record[6]
            start = int(record[3])
            end = int(record[4])
            gene_id = record[8].split(";")[0][3:].replace("TU","model")
            
            # if the orientation of next gene is different from preivous, add last block to the dictionary
            if scaffold_id!=last_scaffold or orientation!=last_orientation:
                if len(new_block) != 0:
                    if new_block.scaffold in blocks:
                        blocks[new_block.scaffold].append(new_block)
                    else:
                        blocks[new_block.scaffold] = [new_block]
                # create a new blcok
                new_block = UniDirectionGeneBlock()
            
            # add the current gene to the block
            new_block.add_gene(gene_id,scaffold_id,start,end,orientation)
            
            last_orientation = orientation
            last_scaffold=scaffold_id

        # this is the end of the for loop, add the last block to the dictionary
        if len(new_block) != 0:
            if new_block.scaffold in blocks:
                blocks[new_block.scaffold].append(new_block)
            else:
                blocks[new_block.scaffold] = [new_block]    
    return blocks


def calculate_gc(seq,baseline=None):
    """
    calcuate the GC% for a given DNA sequence with possible Ns.
    If baseline GC% (usually the background GC%) is provided, it will be used to impute Ns.
    else Ns are igored in GC% calcuation.
    """
    if not baseline:
        count_g = seq.count("G")
        count_c = seq.count("C")
        count_n = seq.count("N")
        if count_n == len(seq):
            return None
        else:
            return (count_g+count_c)/(len(seq)-count_n)
    else:
        count_g = seq.count("G")
        count_c = seq.count("C")
        count_n = seq.count("N")
        return (count_g+count_c+count_n*baseline)/len(seq)

def scan_gc_drop(genome,scaffold,start,end,extension_length,window_size,gc_change_threshold,drop_width_threshold,directory):
    """
    In the given genome interval, we aim to identify a specific region that exhibits a sharp drop in GC% compared to the 
    overall GC% of the interval. i.e. identify "peaks" that represent a significant decrease in GC% compared to the background.

    A sliding window of a specific size was moved across the genome interval, and the GC% was calculated for each window.
    The background GC% was determined as the average GC% of the entire interval. The sliding window analysis identified regions
    where the GC% within the window dropped below the background GC%. This marked the beginning of the "peak."
    The "peak" continued until the GC% within the window increased above the background GC%, indicating the end of the sharp 
    drop in GC%.
    """
    scaffold_seq = genome[scaffold]
    scaffold_len = len(scaffold_seq)
    largest_peak = None
    found = False
    # igore short scaffolds
    if scaffold_len >= extension_length:

        # extend the putative TAD regions by extension length
        scan_start = max(0,start - extension_length//2)
        scan_end = min (scaffold_len,end+extension_length//2)

        peaks = []
        GC_values = []
        # normal_gc stands for GC%>= baseline
        last_subwindow_is_normal_gc = True

        window_seq = scaffold_seq[scan_start:scan_end].seq
        average_gc = calculate_gc(window_seq)
        if not average_gc:
            return None
        baseline_gc = average_gc
  
        for i in range(0,len(window_seq)-window_size):
            subwindow = window_seq[i:i+window_size]
            subwindow_gc = calculate_gc(subwindow,baseline_gc)

            GC_values.append(subwindow_gc)
            if subwindow_gc <= average_gc:
                current_subwindow_is_normal_gc = False
            else:
                current_subwindow_is_normal_gc = True

            # current subwindow gc lower than average
            if not current_subwindow_is_normal_gc:

                # create a peak 
                if  last_subwindow_is_normal_gc:
                    peaks.append([i,1,subwindow_gc,i])
                else:
                    peaks[-1][1] += 1
                    if subwindow_gc < peaks[-1][2]: 
                        peaks[-1][2] = subwindow_gc
                        peaks[-1][3] = i
                    
            last_subwindow_is_normal_gc = current_subwindow_is_normal_gc
            baseline_gc=subwindow_gc
            
    for peak in peaks:
        if peak[1] >= drop_width_threshold and peak[2] <=average_gc-gc_change_threshold:

            peak_start = peak[0] +scan_start
            peak_centre = peak[3] +scan_start

            if peak_centre +window_size//2>= start and peak_centre+window_size//2<= end:
                if not largest_peak:
                    largest_peak = (average_gc-peak[2],peak[1],peak[0],peak_start,GC_values,average_gc)
                elif average_gc-peak[2] > largest_peak[0]:
                    largest_peak = (average_gc-peak[2],peak[1],peak[0],peak_start,GC_values,average_gc)

                found = True
                #return [peak_start, peak[1]]
    

    plt.plot(range(scan_start+window_size//2,scan_start+window_size//2+len(GC_values)),GC_values)
    plt.axhline(average_gc,color='red',ls='--')
    if found:
        hight,peak_len,window_pos,peak_start,GC_values,average_gc = largest_peak       
        plt.axvline(peak_start+window_size//2,color='green',ls='--')
        plt.axvline(peak_start+peak_len+window_size//2,color='green',ls='--')
    plt.axvline(start,color='black',ls='--')
    plt.axvline(end,color='black',ls='--')

    if found:
        plt.title("peak width: %d;peak_hight: %f  "%(peak_len,hight))
        plt.savefig("%s/%s_%d.png"%(directory,scaffold,peak_start))
        plt.clf()
        return True
    else:
        plt.savefig("%s/None_%s.png"%(directory,scaffold))
        plt.clf()
        return None


def main():
    gene_gff_file = sys.argv[1]
    genome_file =  sys.argv[2]
    length_threshold = int(sys.argv[3])  # mininum number of genes in a UniDirectionGeneBlock 
    report_GC_drop = True

    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    blocks = search_unidirectional_gene_blocks(gene_gff_file)

    if not os.path.exists("boundary"):
        os.makedirs("boundary")
    if not os.path.exists("centre"):
        os.makedirs("centre")

    for scaffold in blocks:
        if len(blocks[scaffold]) >= 2:
            for i in range(len(blocks[scaffold])-1):
                current_block = blocks[scaffold][i]
                next_block = blocks[scaffold][i+1]

                out_dir=None
                if len(next_block) >= length_threshold and len(current_block) >= length_threshold:
                    if current_block.orientation == "+":
                        print("boundary found!", end =" ")
                        out_dir="boundary"
                    else:
                        print("centre found!",end=" ")
                        out_dir="centre"

                    start = current_block.end
                    end = next_block.start

                    if report_GC_drop:
                        # scan_gc_drop(genome,scaffold,start,end,extension_length,window_size,gc_change_threshold,drop_width_threshold,directory):
                        if scan_gc_drop(genome,scaffold,start,end,20000,5000,0.05,4000,out_dir):  
                            print("found GC drop!")
                        else:
                            print(" GC drop NOT found!")
                    print(current_block.genes,current_block.start,current_block.end,current_block.orientation)
                    print(next_block.genes,next_block.start,next_block.end,next_block.orientation)
                    print("---------------------")
                    
                
if __name__=="__main__":
    main()