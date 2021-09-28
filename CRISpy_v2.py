#!/usr/bin/env python

from __future__ import division
import os
import glob
import sys
import csv
from collections import Counter
from collections import OrderedDict
import pandas as pd


# CRIS.py
#
# @SeanNesdoly: Version 2 (CRISpy_v2.py) targets Python 2.7; I modified this
# script to work with Python 3 (based on CRISpy_v1-py3.py).

# Modify parameters below in the get_parameters() section.

# [2021-09-09] @SeanNesdoly
# Forked from GitHub repo: https://github.com/patrickc01/CRIS.py
#
# Citation:
#   Connelly J., Pruett-Miller S. CRIS.py: A Versatile and High-throughput
#   Analysis Program for CRISPR-based Genome Editing. Scientific Reports 9,
#   4194 (2019)
#
# patrickc01/CRIS.py is licensed under the 'GNU General Public License v3.0'.
# Link to license: https://github.com/patrickc01/CRIS.py/blob/master/LICENSE

def get_parameters(target_gene):
    # Note to user: change text inside of quote marks ('YOUR DNA SEQUENCES GO
    # HERE') for your experiment. Case of text does not matter.

    # Name of output directory (see 'save_dir'); used for calculating
    # sgRNA-specifc 'raw_wt_counter' values
    ID = target_gene

    if (ID == 'DDX6'):
        ref_seq   = 'TCATTAAGCAGCTCAGGACTGTAATATTATAAACCTTTATCCTCTTTGCTCAAATTAAAATTAATGAGAATGTATGTTTTCTAAATCTTGCATAGTGCACGTGGTGATTGCTACCCCTGGGAGAATCCTGGATCTTATTAAGAAAGGAGTAGCAAAGGTTGATCATGTCCAGATGATAGTATTGGATGAGGTAATGTCTCTTCATTTGGCATTATACTGTTCTGTACTTTTCTGTGCTTTGCTTATAACTGCCAGTGAATGAAGCAAAAGCTGCTTGTGTTTATGATCTAATAAACACATAAACAATCAGTTTAAATTAATTTGCTTATAAAATGTGTTTTACCCGTCAGAGGG'
        seq_start = 'CATTAAGCAGCTCAG'
        seq_end   = 'ATGAGGTAATGTCTC'
    elif (ID == 'SMG9'):
        ref_seq   = 'CAGCTCACAGGGTCACACTTCGGATGCATCCAAGGCTGAGAACCAGAGCTCTCACACACCCCTGCCTGCTCCTCTCCACAGGACTTATGTTTTCCGGGCCCAGAGCGCTGAAATGAAGGAACGAGGGGGCAACCAGACCAGTGGCATCGACTTCTTTATTACCCAAGAACGGATTGTTTTCCTGGACACACAGGTGCCAGCCCCGCCTGCCTTTGCCCCACCCTGCCACATCCTGCATTGAGCCTCTTCCTGGAATGTGATTAAGCCAGAACCCTGCGTGGATTCCATGTGCTTCCAGTCAAGCCAT'
        seq_start = 'AGCTCACAGGGTCAC'
        seq_end   = 'ACACACAGGTGCCAG'
    elif (ID == 'CARM1'):
        ref_seq   = 'AGGAGTGCAGGAACGAATGGATGACAGGCTGGGAGCACCCAGGGTTGGGGGTCTTGGGGTCCTTTGGAAGCTCTGCCAGGCAGTGTGGGATGTGTCACCTGACGCCAGCACCCCTCCCTGCCCCACTCCCAGGAAACATGTTTCCTACCATTGGTGACGTCCACCTTGCACCCTTCACGGATGAACAGCTCTACATGGAGCAGTTCACCAAGGCCAACTTCTGGTGAGTGTGCCCTGGGTGTCCCGCCTGGGCCCCACAGCCTGCCTTCTCAGGGACAGCCCCAGCTCCCCAGAGAGCCTGCACTCCTCTTTTTCTGAAAGACTTGGGCTAGATGAGGGC'
        seq_start = 'GGAGTGCAGGAACGA'
        seq_end   = 'CCAAGGCCAACTTCT'
    elif (ID == 'NXF1'):
        # As the reverse complemented sgRNA sequence for NXF1 was found in the
        # reverse reads (*R2*.fastq.trimmed), the below sequences have been
        # reverse complemented to match this.
        ref_seq   = 'TGGGAGCTACTGGGTCCTTGTATGGCTGTAGACCACAGAGCCCATACTTTCTTTTTTGTTTTTCCCATAGGTACTATGCAATTTACGACTCTGGAGACCGACAAGGGCTCCTGGATGCCTACCATGATGGGGCCTGCTGTTCCCTGAGCATTCCTTTCATTCCTCAGAACCCTGCCCGGTTAGTATCACATCCCAGAATCTGCCCGGGACCACTGGGTTTCCCAGCAGATACAGGGCACACCTGCTGGCAGCACTGCTCTGAACACCACGTCCATTCCTTTTATTCCCCCAGAAGCAGCTTAGCCGAGTA'
        seq_start = 'GGGAGCTACTGGGTC'
        seq_end   = 'TCACATCCCAGAATC'
    elif (ID == 'NXT1'):
        ref_seq   = 'TGTTCCATTGTCAGGGCAGAATGAACTTTGGCATTCACGTGGCTTCTCTTCAACCTTACTTCCCTGCAGCCCCTGGTTCCCCAAGGCAGAGGAAATACCCTGGTGGAGCCCTCCTTCCATAGAACCAGAGATGGCATCTGTGGATTTCAAGACCTATGTGGATCAGGCCTGCAGAGCTGCTGAGGAGTTTGTCAATGTCTACTACACCACCATGGATAAGCGGCGGCGTTTGCTGTCCCGCCTGTACATGGGCACAGCCACCCTGGTCTGGAATGGCAATGCTGTTTCAGGACAAGAATCCTTGAGTGAGTTTTTTGAAATGTTGCCTTCCAGCGAGTTCCAAATCAGCGTGGTAGACTGCCAGCCTGTTCATGATGAAGCCACACCAAGCCAGACCACGGTCCTTGTTGTCATCTGTGGATCAGTGAAGTTTGAGGGGAACAAACAACGGGACTTCAACCAGAACT'
        seq_start = 'GTTCCATTGTCAGGG'
        seq_end   = 'CTGGTCTGGAATGGC'

    test_list = [
        str('DDX6'), str.upper('ACGTGGTGATTGCTACCCCT'),
        str('SMG9'), str.upper('GCTGAAATGAAGGAACGAGG'),
        str('CARM1'), str.upper('TTCACGGATGAACAGCTCTA'),
        str('NXF1'), str.upper('ACGACTCTGGAGACCGACAA'),  # reverse complement
        str('NXT1'), str.upper('ACTACACCACCATGGATAAG')
    ]

    # Choose R1 or R2 reads for analysis based on target gene (ID)
    if (ID == 'NXF1'):
        fastq_files = '*R2*.fastq.trimmed' # use reverse reads
    else:
        # Use forward reads for genes 'DDX6', 'SMG9', 'CARM1', & 'NXT1'
        fastq_files = '*R1*.fastq.trimmed'

    # Command-line argument parsing, added by ariva@ufl.edu
    cmdline_fastqs = []
    args = sys.argv[1:]
    if "-h" in args or "--help" in args:
        usage()
    prev = ""
    for a in args:
        if prev in ["-i", "--id"]:
            ID = a
            prev = ""
        elif prev in ["-r", "--ref"]:
            ref_seq = str.upper(a)
            prev = ""
        elif prev in ["-s", "--start"]:
            seq_start = str.upper(a)
            prev = ""
        elif prev in ["-e", "--end"]:
            seq_end = str.upper(a)
            prev = ""
        elif prev in ["-t", "--test"]:
            test_list = read_test_list(a)
            prev = ""
        elif a in ["-i", "--id", "-r", "--ref", "-s", "--start", "-e", "--end", "-t", "--test"]:
            prev = a
        else:
            cmdline_fastqs.append(a)

    if cmdline_fastqs:
        fastq_files = cmdline_fastqs
    else:
        data_dir = os.getenv("HOME") \
            + '/scratch/210721_MultiKO_CellLineSeqValidation/out/trimmed_reads/06_m50bp/rcFullSeq_rc19bpSuffix_SamplesRemoved/'
        fastq_files = glob.glob(data_dir + fastq_files)

    # End additions
    print("Input files:")
    for f in fastq_files:
        print("  " + f)
    return ID, ref_seq, seq_start, seq_end, fastq_files, test_list

# Added by ariva@ufl.edu
def read_test_list(filename):
    """Read the test_list from a tab-delimited file `filename'. The first
column contains the test name, the second column contains the test sequence."""
    tl = []
    with open(filename, "r") as f:
        c = csv.reader(f, delimiter='\t')
        for line in c:
            tl.append(line[0])
            tl.append(str.upper(line[1]))
    return tl

def usage():
    sys.stdout.write("""CRIS.py [options] fastqs...
Where options are (short form followed by long form):
  -h   | --help    | Display this help message.
  -i I | --id I    | Set run ID to I.
  -r R | --ref R   | Set reference sequence to R.
  -s S | --start S | Set start sequence to S.
  -e E | --end E   | Set end sequence to E.
  -t T | --test T  | Read test list from tab-delimited file T.
The tab-delimited file passed as the value of `-t' should have two columns,
containing the test name and test sequence respectively.
""")
    sys.exit(0)

def pairwise(iterable):
    #Make an ordered dictionary from items in in test list
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a) # @SeanNesdoly: izip from itertools is zip in Python 3


def make_project_directory(save_dir):
    #Create a project directory to store files in
    cwd = os.getcwd()   #Get current working directory
    try:
        os.stat(save_dir)
    except:
        os.mkdir(save_dir)

def write_to_file(record_entry,f):
    #Writes results to the results_counter.txt file. Method of inserting space
    #between each well to make the .txt file easier to read. Based on # of
    #items in well list
    f.write("\n")
    for record in record_entry:         #record entry is a giant list (master_Record) that contains well-lists as elements
        for line in record:
            if len(record) > 2:
                f.write(str(line) + '\n')
            else:
                pass
        if len(line) > 2:
            f.write('\n\n')
        else:
            pass

def make_counter(indel_size_list, current_fastq_file, dict_Counters, c_Counter, SNP_test, raw_wt_counter):
    top_common = 12             #Number of top found reads to write to results_counter .txt file
    temp_counter = Counter(indel_size_list).most_common(top_common)    #Count top indels present in fastq file
    temp_dict = OrderedDict()

    #If there are not a total of at least 'top_common' common reads from the
    #summary file, fill the spaces with NA, NA
    if len(Counter(indel_size_list).most_common(top_common)) < top_common:
        for i in range (0, top_common - len(Counter(indel_size_list).most_common(top_common))):
            temp_counter.append(("NA", "NA"))
    else:
        pass

    temp_dict['Name'] = str(current_fastq_file)              #Fill in Name column of CSV file with the file fastq name
    temp_dict['Sample'] = ''                                 #Makes a blank column for user to fill in sample names

    for k in dict_Counters:
        try:
            temp_dict[k] = str(str(dict_Counters[k]) + ' (' + str(format((dict_Counters[k] / c_Counter*100), '.1f')) + '%)')
            temp_dict[str('%'+ k)] = str(format((dict_Counters[k] / c_Counter)*100, '.1f'))
            temp_dict['Total'] = c_Counter
            temp_dict['SNP_test'] = SNP_test
            temp_dict['raw_wt_counter'] = raw_wt_counter
            temp_dict['Total_indel'] = str((format(((c_Counter - indel_size_list.count(0))/c_Counter)*100,'.1f')))
        except ZeroDivisionError:
            pass

    counter = 1
    for k,g in temp_counter:        #Fill in top indels sizes and amount
        temp_dict['#'+ str(counter)+ '-Indel']=k
        try:
          # temp_dict['z%Indel'+str(counter)]= format(g/c_Counter *100, '.1f')
           temp_dict['#'+ str(counter)+'-Reads(%)'] = str(g) + ' ('+ str(format(g/c_Counter*100, '.1f')) + '%)'
        except TypeError:
            pass

        counter+=1

    return temp_dict



def search_fastq(ID, ref_seq, seq_start, seq_end, fastq_files, test_list):
    #Process the fastq file and look for reads
    test_dict = OrderedDict()     #Name, Sequence of each item that is being searched for
    dict_Counters = OrderedDict() #Name, Count of each item in the dictionary being searched for (ie current_fastq_file)
    fastq_counter = 0             #Count the number of fastq_files with reads
    master_Record = []            #Master list, Master list contains lists
    indel_size_list = []          #list of indel sizes found in the fastq.  For each read, one indel size is recorded
    csv_summary_df = pd.DataFrame()
    master_distance_and_count_summary = []

    for x,y in pairwise(test_list):          #Create an ordered dictionary of items in test_list
        test_dict[x] = y

    save_dir = os.getenv("HOME") \
        + '/scratch/210721_MultiKO_CellLineSeqValidation/out/crispy/' \
        + str(ID) \
        + '/'
    make_project_directory(save_dir)
    print("Working directory: {}".format(str(os.getcwd())))
    print("Output directory: {}".format(save_dir))

    file_name = save_dir + 'results_counter_' + ID + '.txt'
    with open(file_name, "w") as f:
        wt_distance = ref_seq.find(seq_end)+len(seq_end) - ref_seq.find(seq_start)      #Expected size of WT read, the distance between the two anchor points made from seq_start and seq_end
        f.write("Starting CRISpy analysis for: " + ID + '\n')
        print("Starting CRISpy analysis for: " + ID + '\n')
        f.write(str("seq_start: "+seq_start+'\n'))
        f.write(str("seq_end: "+seq_end+'\n'))
        f.write("Test_Sequences (from 'test_list'): \n")
        for key, value in test_dict.items():                #Go through the test_dict and write each item that is being searched for
            f.write(str(key)+": "+value+'\n')
            print(key, value)
            dict_Counters[str(key)]=0

        print('Expected WT distance: {}'.format(wt_distance))     #The expected distance between  seq_start and seq_end if the DNA is WT/ REF
        if wt_distance < 0:
            f.write("\n\n WARNING: THIS IS NOT GOING TO GIVE YOU THE FULL DATA. YOUR EXPECTED WT DISTANCE IS LESS THAN 0, it is: {}\n  Check your seq_start and seq_end again\n".format(wt_distance))
        else:
            pass

        print("Program Running")

        for each_fastq_file in fastq_files:   #For each clone (really each fastq file in directory), open the file as "clone"
            c_Counter = 0                     #Reset control counter to 0, this counter counts how many times both seq_start and seq_end are found in a line.
            start_counter = 0                 #How many times seq_start is found in a fastq file, used for SNP check
            end_counter = 0                   #How many times seq_end is found in a fastq file, used for SNP check
            raw_wt_counter = 0                #Calculate how many times first item in test_list is found in fastqfile, a check for SNPs

            for items in test_dict:            # RESET DICT COUNTERS TO 0
                dict_Counters[str(items)] = 0

            indel_size_list = []                          #List of all the indel sizes found in each file
            fastq_name = str(each_fastq_file)             #Well name that contains the clone being screened
            line_list =[]                                 #List of all the lines found in a fastq file.  These lines must contain both seq_start and seq_end.  Used to find most highly occuring sequences

            with open(str(each_fastq_file), "r") as current_fastq_file:
                for line in current_fastq_file:   #For each line in the fastq file
                    #Counts the number of times the first item in test_dict
                    #is found in ANY of the lines of the fastq file. This
                    #is a check in case SNPs are in BOTH seq_start and
                    #seq_end

                    # @SeanNesdoly: use the sequence of the sgRNA currently
                    # under analysis, instead of arbitrarily using the
                    # first one in 'test_dict'. This will give more accurate
                    # sgRNA-specific counts.
                    if test_dict[ID] in line:
                        raw_wt_counter += 1
                    if line.find(seq_start)>0 and line.find(seq_end)>0:
                        c_Counter += 1
                        start_counter += 1 #Count # of times seq_start is found
                        end_counter += 1   #Count # of times seq_end is found
                        read_start = line.find(seq_start)
                        read_end = line.find(seq_end)+len(seq_end)
                        indel_size = line.find(seq_end)+len(seq_end) - line.find(seq_start) - wt_distance
                        indel_size_list.append(indel_size)
                        line_list.append(line[read_start:(read_end)])
                        for item in test_dict:
                            if test_dict[item] in line:
                                dict_Counters[item] += 1
                            else:
                                pass
                    elif line.find(seq_start)>0 and line.find(seq_end)<0:        #If seq_start is found and seq_end is not found, for SNP test
                        start_counter += 1
                    elif line.find(seq_end)>0 and line.find(seq_start)<0:        #If seq_end is found and seq_start is not found, for SNP test
                        end_counter += 1
                    else:
                        pass

            try:
                SNP_test = format(start_counter / end_counter, '.2f')        #Compare counts of seq_start and seq_end for SNP test
            except ZeroDivisionError:
                pass
            try:
                #Calculate raw_wt_counter
                raw_wt_counter = str(str(raw_wt_counter) +  ' ('+ format((raw_wt_counter/ dict_Counters[ID]), '.1f') +')')
            except ZeroDivisionError:
                pass

            if c_Counter == 0:
                pass
            elif c_Counter > 10:  #if more than 10 control read counts, record data
                print("{}: Total_reads:{}, {}".format(fastq_name, str(c_Counter).ljust(2), dict_Counters.items()))
                fastq_counter += 1

                test_list_string=str(" Testing: ")
                for k,v in dict_Counters.items():
                    test_list_string=test_list_string+"({}:{}), ".format(k,v)

                temp = Counter(line_list).most_common(12)
                #summary_line is a list, format:  Miller-Plate13-C01 TOTAL:2072 OrderedDict([('g3', 2010), ('Block_Only', 0), ('Mod_Only', 2), ('Block_Mod', 0), ('Full', 0)])       [(0, 2070), (-1, 2)]
                summary_line = ([str(fastq_name) + " TOTAL:" + str(c_Counter)+" "+test_list_string+"     "+"Top_reads:"+ str(Counter(indel_size_list).most_common(12))])

                for k,v in temp:          #Append the top found DNA sequences to summary_line
                    summary_line.append('{} , {}'.format(k,v))

                master_Record.append(summary_line)
                master_distance_and_count_summary.append(make_counter(indel_size_list,
                                                                      str(fastq_name),
                                                                      dict_Counters,
                                                                      c_Counter,
                                                                      SNP_test,
                                                                      raw_wt_counter))
            else: # c_Counter > 0 && c_Counter <= 10 (cannot be < 0)
                 pass

        print("SUMMARY")
        make_project_directory(ID)
        #print(master_distance_and_count_summary)

        pd_columns = ['Name', 'Sample', 'Total', 'Total_indel', '#1-Indel',
                      '#1-Reads(%)', '#2-Indel', '#2-Reads(%)', '#3-Indel',
                      '#3-Reads(%)', '#4-Indel', '#4-Reads(%)', '#5-Indel',
                      '#5-Reads(%)', '#6-Indel', '#6-Reads(%)', '#7-Indel',
                      '#7-Reads(%)', '#8-Indel', '#8-Reads(%)', 'SNP_test',
                      'raw_wt_counter']

        flip_dict = list(test_dict.items())   #Need to flip order of items in dictionary.  That way when they are inserted into the excel list, the order will come out correct
        flip_dict.reverse()
        flip_dict = OrderedDict(flip_dict)

        for k, v in flip_dict.items():        #Insert the items from test_dict into position 3 for the column output, after 'Total'
            pd_columns.insert(3, k)

        for k, v in test_dict.items():        #Insert the % values for test_dic at the end
            pd_columns.append(str('%'+k))

        csv_summary_df = pd.DataFrame.from_records(master_distance_and_count_summary, index='Name', columns=pd_columns)
        csv_summary_df = csv_summary_df.sort_index()
        csv_summary_df = csv_summary_df[pd.notnull(csv_summary_df['Total'])]

        try:
            csv_summary_df.to_csv(str(save_dir+ID)+'.csv')    #Filename to save csv as
        except (IOError):
            print('ERROR.  Script did not execute properly.')
            print(('The requested .csv file {} is either open or you do not have access to it.  If open, please close file and rerun program').format(str(save_dir+ID)+'.csv'))

        master_Record = sorted(master_Record)
        print("Total wells with product:", fastq_counter)
        write_to_file(master_Record, f)

def main():
    ID = ''
    ref_seq = ''
    seq_start = ''
    seq_end = ''
    fastq_files = ''
    test_list = []

    print("CRISpy_v2.py\nModified by @SeanNesdoly\nMain program")

    # @SeanNesdoly: Run CRISpy analysis on each target sgRNA/gene in turn.
    genes = ['DDX6', 'SMG9', 'CARM1', 'NXF1', 'NXT1']
    for g in genes:
        ID, ref_seq, seq_start, seq_end, fastq_files, test_list = get_parameters(g)
        search_fastq(ID, ref_seq, seq_start, seq_end, fastq_files, test_list)
        print("Completed CRISpy analysis for: " + g + '\n\n')

    print("Done")

if __name__== "__main__":
    main()
