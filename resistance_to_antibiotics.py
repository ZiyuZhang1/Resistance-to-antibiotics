#!/usr/bin/env python3
import sys, re, gzip,time

time_start = time.time()

### functions ###
# global data declear
reference_kmer_set = set()
reference_kmer_set_dict = dict()
reference_kmer_dict = dict()

def cut_kmer(seq):
    '''cut squence into kmers(K=19), and create a list of kmers'''
    kmer_list = list()
    for i in range(len(seq)-18):
        kmer = seq[i:i+19]
        kmer_list.append(kmer)
    return kmer_list


def read_reference_file(reference_filename):
    '''Function to read reference file and create three databases for reference gene(reference_kmer_set(contains all reference gene kmers), reference_kmer_set_dict(ID as key, kmer set as value) and reference_kmer_dict{ID:[[kmer, kmer, ...],[[0]*len(seq)]})'''
    try:
        infile = open(reference_filename, 'r')
        # find first header line
        line = infile.readline()
        while line != '' and line[0] != '>': line = infile.readline()
        # read the rest of the fasta file
        while line != '':
            # get the header as ID
            ID = line[:-2]
            # read sequence until next header or end of file
            seq = ''
            line = infile.readline()
            while line != '' and line[0] != '>':
                line = re.sub("\s", "", line)
                seq += line
                line = infile.readline()
            # after the sequence is read, now create three databases for reference gene
            if ID in reference_kmer_dict:
                print("dulplicated gene: ", ID)
            else:
                # create reference_kmer_dict{ID:[[kmer, kmer, ...],[[0]*len(seq)]}
                kmer_list = cut_kmer(seq)
                kmer_infor = list()
                kmer_infor.append(kmer_list)
                kmer_infor.append([0]*len(seq))
                reference_kmer_dict[ID] = kmer_infor
                # create reference_kmer_set(contains all kmers)
                kmer_set = set(kmer_list)
                reference_kmer_set.update(kmer_set)
                # create reference_kmer_set_dict(ID as key, kmer set as value)
                reference_kmer_set_dict[ID] = kmer_set   
        infile.close()
    except IOError as error:
        sys.stderr.write("File I/O error, reason: " + str(error) + "\n")
        sys.exit(1)
    except ValueError as error:
        sys.stderr.write("Format error: " + str(error) + "\n")
        sys.exit(1)

def map_reads(sample_filename):
    '''map reads to resistance gene'''
    try:
        infile = gzip.open(sample_filename, 'rt')
        # find reads line
        line = infile.readline()
        seqflag = False
        line_number = 0
        while line != '':
            if line[0] == '+':
                seqflag = False
            if seqflag:
                line = re.sub("\s", "", line)
                read = line
                map(read)
            if line[0] == '@':
                seqflag = True
            line = infile.readline()
            line_number += 1
        infile.close()
    except IOError as error:
        sys.stderr.write("File I/O error, reason: " + str(error) + "\n")
        sys.exit(1)
    except ValueError as error:
        sys.stderr.write("Format error: " + str(error) + "\n")
        sys.exit(1)
    
def map(read):
    '''map read to reference gene dictionary, if all mapped, add one to every kmer count'''
    # filter and precise match
    strand_flag = False
    # condition 1: read matched reference strand
    if not strand_flag:
        # roughly filter reads whose most kmers don't contain in reference_kmer_set
        # if read_set is not different from referemce_kmer_set, means it is possible they can match
        read_list = cut_kmer(read)
        read_set = set(read_list) 
        if not read_set.isdisjoint(reference_kmer_set):
            # determine whether read can match with this reference gene
            for ID in reference_kmer_set_dict:
                if not read_set.isdisjoint(reference_kmer_set_dict[ID]):
                    matched_ID = ID
                    result = precise_map(read,matched_ID)
                    if result is not None:
                        for i in range(result[0],result[1]):
                            reference_kmer_dict[matched_ID][1][i] += 1
                    strand_flag = True
    # consider complement strand
    else:
        # translate read
        complementtable = str.maketrans('ATCGatcg', 'TAGCTAGC')
        read = read.translate(complementtable)[::-1]
        # roughly filter reads whose most kmers don't contain in reference_kmer_set
        # if read_set is not different from referemce_kmer_set, means it is possible they can match
        read_list = cut_kmer(read)
        read_set = set(read_list) 
        if not read_set.isdisjoint(reference_kmer_set):
            # determine whether read can match with this reference gene
            for ID in reference_kmer_set_dict:
                if not read_set.isdisjoint(reference_kmer_set_dict[ID]):
                    matched_ID = ID
                    result = precise_map(read,matched_ID)
                    if result is not None:
                        for i in range(result[0],result[1]):
                            reference_kmer_dict[matched_ID][1][i] += 1

def precise_map(read,matched_ID):
    '''precise match this read to reference_kmer_dict[matched_ID], change the count of reference_kmer_dict[matched_ID]'''
    # if try head match, end match is meanless
    try_end_match = False
    for index, kmer in enumerate(reference_kmer_dict[matched_ID][0]):
        # initialise variables
        is_match = False
        start = 0
        end = 0
        if read[0:19] == kmer:
            is_match = True
            try_end_match = False
            start = index
            end = index + len(read)
            # check if all base matched
            for i in range(19, len(read)-18):
                try: 
                    if read[i:i+19] != reference_kmer_dict[matched_ID][0][index+i]:
                        is_match = False
                        try_end_match = True
                        break
                # only end parts of reads doesn't matched
                except IndexError as error:
                    end = len(reference_kmer_dict[matched_ID][0]) + 18
                    try_end_match = False 
                    break 
       
        if try_end_match:
            if read[-19:] == kmer:
                is_match = True
                end = index + 19
                start = 0
                # check if all base matched
                for i in range(len(read)-end,len(read)-18):
                    if read[i:i+19] != reference_kmer_dict[matched_ID][0][i-(len(read)-end)]:
                        if index-i < 0:
                            is_match = False
                            break

        if is_match:
            return [start,end]
            
### Main code ###
# get the file
if len(sys.argv) == 1:
    reference_filename = input('please enter a reference filename: ')
    sample_filename = input('please enter a sample gzipped filename: ')
elif len(sys.argv) == 3:
    reference_filename = sys.argv[1]
    sample_filename = sys.argv[2]
else: 
    sys.stderr.write('Usage: complement.py<filename>\n')
    sys.exit(1)

# use lists of dicts to make two database of resistance gene, mapped_dict and reference_kmer_dict
read_reference_file(reference_filename)
# map reads to resistance gene
map_reads(sample_filename)

# open outfile
outfile = open("mapped_result_file", "w")
# calculate coverage and depth; write into outfile
for ID in reference_kmer_dict:
    depth_sum = 0
    matched_kmer_count = 0
    depth_sum = sum(reference_kmer_dict[ID][1])
    if depth_sum > 0:
        for base_count in reference_kmer_dict[ID][1]: 
            if base_count > 0:
                matched_kmer_count +=1
        coverage = matched_kmer_count/len(reference_kmer_dict[ID][1])
        depth = depth_sum/matched_kmer_count
        if coverage > 0.95 and depth > 10:
                outfile.write(ID + "\t"+ "coverage"+ str(coverage) + "\t"+ "depth" + str(depth)+ "\n")
        # outfile.write(ID + "\t"+ "coverage"+ str(coverage) + "\t"+ "depth" + str(depth)+ "\n")
outfile.close()
        # print(ID, "\t", "coverage", str(coverage), "\t" + "depth", str(depth) + "\t")
        # if coverage > 0.95 and depth > 10:
        #     outfile.write(ID + "\t"+ "coverage"+ str(coverage) + "\t"+ "depth" + str(depth)+ "\t")

time_end = time.time()
print("time cost", time_end-time_start, "s")