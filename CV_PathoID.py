#!/usr/bin/python

"""##################################################################
# Script to automatically find and retrieve a list of variants
#   on ClinVar, for their clinical significance and (if present)
#   associated disease condition.
#
# Harcoded reading input file format:
#   - column 10 (index 9): Gene symbol
#   - column 8 (index 7): Variant function type (VFT)
#   - column 11 (index 10): detailed variant annotation
#   - column 12 (index 11): All SNPs (e.g. rs number)
#   * For future: maybe attempt to recognize header string in order
#       to automatically detect which columns are of interest
#   * Also, be careful that the interest_cols list order and the
#       Var() constructor input order matches
#
#
# Note to self:
#   - sys.argv[1] is the input file path
#   - Only identifies and reads .csv and .output files for now
#   - The interestd column numbers are hard-coded in read_file
#   - Problem with newline - not marked by '\n', will need to
#       manually separate lines via delimiters "\r"
#       - This is assumed to be the same for both .output and .csv
#           files, which may or may not cause problems down the road
#       - However, for now, this seems to be the case and there has not
#           been any problems running either .csv or .output files yet
#   - Outputting an appended file:
#       - Re-read and split the file, then re-write with the pathogenic
#           and disease condition status appended to end of each line
#       - Should probably error check for match
#
# Potential ideas:
#   - Implement pickle files in order to separate longer search processes
#       - i.e. check for pre-existing pickle files, separate each
#           search down to small steps to avoid network errors and such
#
#
# Author: Anthony Chen
##################################################################"""
#System libraries
import sys
import re
import datetime
import csv
#User-created files
from variant import Var
import variant
import connect

"""### Functions to open the correct type of file and initialize ###"""
#Function that oversees file input
def read_file(variant_list, filename, interest_cols):
    #Initiate the array containing Str data from the file
    file_content=[]
    #Identify the type of spreadsheet file via regex
    if re.search(r'\.output$', filename):
        #Initialize tab-delimited file format
        print "[Program] Initializing tab-delimited .output file..."
        if (read_and_initialize(filename, '\t', interest_cols, file_content) != 0):
            return 1
    elif re.search(r'\.csv$', filename):
        #Initialize comma-separated file format
        print "[Program] Initializing comma-separated .csv file..."
        if (read_and_initialize(filename, ',', interest_cols, file_content) != 0):
            return 1
    else:
        print "[ERROR] File extension not recognized, make sure it is 'output' or 'csv'"
    #Initialize the content of the file:
    for r in file_content:
        variant_list.append(Var(r[0],r[1],r[2],r[3]))
    #Return good
    return 0

# Function that opens and reads from a file, then initialize and returns
# a 2D list containing the pre-specified columns of interest
def read_and_initialize (filename, delim, interest_cols, file_content):
    #Attempt to open the input stream
    try:
        f = open(filename, 'r')
    #Except for I/O error in case file is not present
    except IOError as e:
        print "[ERROR] I/O error({0}): {1}".format(e.errno, e.strerror)
        return 1
    #Read the entire file and separate by newline delimiters
    rows = f.read().split('\r')
    rows.pop(0) #Remove the header (first) row!
    #Iterate through each now-separated row
    for r in rows:
        #Initiate a list containing the current column
        col_content = []
        #Separate each row into columns
        cols = r.split(delim)
        #Iterate through the list of column numbers of interest
        for i in interest_cols:
            col_content.append(cols[i])
        #Append the content into the returnable content array
        file_content.append(col_content)
    #Close the input stream and return good
    f.close()
    return 0


"""################## Function to search ##################"""
def search_ClinVar(variant_list):
    #Edit the annotations to become searchable
    print "[Program] Formatting variant annotations..."
    variant.format_variantList(variant_list)
    #Loop through the variant list to find record IDs via eSearch
    print "[Program] Beginning ClinVar eSearch for record IDs..."
    connect.ClinVar_Search_Loop(variant_list, 0)
    #Loop through the variant list to find record information
    print "[Program] Beginning ClinVar eSummary for record information..."
    if ( connect.ClinVar_Search_Loop(variant_list, 1) != 0):
        return 1
    #If done, return success
    return 0

"""################## Output Functions ##################"""
#Master file to write output
def write_output_file(filename, v_list, output_type):
    #If the user only wants the short .csv summary file
    if output_type == 0:
        print "[Program] Writing .csv ClinVar result summary file..."
        write_new_csvFile(filename, v_list)
    #If the user wants to output a file with the ClinVar info appended
    elif output_type == 1:
        print "[Program] Appending ClinVar result to new copy of original file..."
        append_end_column(filename, v_list)
    #If the user wants both of the above files
    else:
        print "[Program] Writing .csv ClinVar result summary file..."
        write_new_csvFile(filename, v_list)
        print "[Program] Appending ClinVar result to new copy of original file..."
        append_end_column(filename, v_list)
    #Prompt user after finishing file writing and return success
    print "[Program] Done writing output files."
    return 0

#Ask the user for what type of output they want
def get_output_type():
    #Indicate to user file types
    out_types = ["0 - Short .csv file summarizing ClinVar results"]
    out_types.append("1 - Same file as inputted, with appended ClinVar results")
    out_types.append("2 - Both of the above outputs")
    print "[Program] Please indicate the type of output file(s) you wish to create:"
    print "\t%s\n\t%s\n\t%s"%(out_types[0],out_types[1],out_types[2])
    #Prompt the user until a valid user input is received
    while (True):
        usr_in = raw_input('Type [0/1/2] followed by "enter": ')
        #Check if user input is a single digit between 0-2
        if len(usr_in)==1 and usr_in.isdigit() and int(usr_in) < 3:
            #Let the user know the file type they have selected
            print "[Program] You selected: %s" % (out_types[int(usr_in)])
            #Break infinite loop
            break
        #If not valid, prompt user to enter again
        print "[ERROR] Invalid input, try again."
    #Return the integer representation of the user's input
    return int(usr_in)

#Function to generate a new .csv file containing the essential information
def write_new_csvFile(filename, v_list):
    #Generate output file name using the input file name and today's date
    date = str(datetime.date.today()).replace('-','')
    output_name = filename.split('.')[0]+'_ClinVarResultSummary_%s.csv'%(date)
    #Open the output file
    f_out = open(output_name, 'w')
    #Write header
    header = ['Gene Name','Detailed Variant Annotation','SNP (rs number)','Clinical Significance','Disease Conditions']
    f_out.write(','.join(header)+'\r')
    #Iterate through all variants to write content
    for v in v_list:
        #Get all related information
        infoList=[v.gene,v.annotation,v.snp,v.output_clin_sig(),v.output_conditions()]
        #Write information
        f_out.write(','.join(infoList)+'\r')
    #Close the file output stream
    f_out.close()


#Function to append the pathogenicity results to the last column
# *Note, the file type is preserved
def append_end_column(filename,v_list):
    #Set the delimiter
    if re.search(r'\.output$', filename): #tab-delimited .output file
        delim = "\t"
    elif re.search(r'\.csv$', filename): #csv file
        delim = ","
    #have another else block?
    #Open reading file - should be okay with exception since it is the same file
    f_in = open(filename, 'r')
    #Read the input file and split via newline delimiter
    input_rows = f_in.read().split('\r')
    #Generate output file name using the input file name and today's date
    date = str(datetime.date.today()).replace('-','')
    output_name = filename.replace('.','_ClinVarAppended_%s.'%(date))
    #Open output file
    f_out = open(output_name, 'w')
    #Write the header of output file
    f_out.write(input_rows[0]+delim+"Clinical Significance"+delim+"Conditions"+'\r')
    #Iterate through the rest of the file
    for i in range(1, len(input_rows)):
        #Find and double check variant information
        if v_list[i-1].annotation not in input_rows[i]:
            print "We have a problem"
            #DO SOMETHING HERE!!!!!!!
        #If good, output the formatted searched information
        patho = v_list[i-1].output_clin_sig()
        cond = v_list[i-1].output_conditions()
        #Write out the pre-existing file content
        f_out.write(input_rows[i])
        #Write out the pathogenicity content
        f_out.write('%s%s%s%s\r'%(delim,patho,delim,cond))
    #Close the input and output files
    f_in.close()
    f_out.close()




# Main Workflow Function
def main():
    #Harcoded columns of interest to initiate from the file
    cols_of_interest = [9,7,10,11]
    #Check for number of command-line arguments
    if len(sys.argv) != 2:
        print "[ERROR] Incorrect number of arguments."
        return 0
    #Ask the user about output options
    output_type = get_output_type()
    #Initiate a list that will contain variant objects
    variant_list = []
    #Open the file and initiate content
    if ( read_file(variant_list, sys.argv[1], cols_of_interest) != 0):
        print "[ERROR] Something went wrong while reading the file"
        return 0
    #Step to make sure that the variant list isn't empty
    print "[Program] %d variants initiated" % (len(variant_list))
    if ( len(variant_list) == 0):
        print "[ERROR] No variants initiated from file"
        return 0
    #Search step
    if ( search_ClinVar(variant_list) != 0 ):
        return 0
    #Output step
    write_output_file(sys.argv[1], variant_list, output_type) #add more error checks?


main()
