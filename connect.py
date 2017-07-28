#!/usr/bin/python

"""##################################################################
# Helper class to connect to ClinVar via the internet
#
# Note to self:
#   - Separating the eSearch and eSummary parts into two separate
#       searches / iterations through all variants
#   - Check for further errors in the eSearch step? (e.g. network error)
#       - Currently, a good type will be returned as long as there
#           is an "IdList" element in the xml
#   - Return exit status from returned results?
#   - Add search indicators (how many searched vs. how many left)
#   - Does not provide the list of conflicting conditions for "conflicting
#       interpretations of pathogenicity" (for eSummary)
#       - either do manual check or future additional steps to add info
#   - Exception handling for clin_sig and conds:
#       - If no text: either will print "None" or attribute error
#       - Not sure if attribute error is the only one that will be thrown
#       - Even print an error message for conditions not found? I think it
#           is a rather common occurance
#       - Likely should do more error handling in the eSummary getReuslt
#
# Author: Anthony Chen
##################################################################"""
import urllib
import time
import xml.etree.ElementTree as ET

#Overall search loop that iterates over all variant records to perform search,
#   while maintaining the appropritate request rate (< 3 / second)
#       Search types:
#           0 - searches for presence of record and retrieve IDs
#           1 - searches for summary and pathogenicity
def ClinVar_Search_Loop(v_list, search_type):
    #Timer variables to prevent excessive ClinVar requests
    start_time = time.time()
    request_counter = 0

    #Iterate through each variant to initiate search terms
    for i in range(0, len(v_list)):
        #User prompt for search progress
        if i % 10 == 0:
            print "\t%d of %d variants"%(i+1,len(v_list))

        #If the current variant is not searchable, skip to the next one
        if v_list[i].searchable == False: continue

        #Generate the query term for the appropriate search
        if search_type == 0: #eSearch query term
            url_query = eSearch_generate_query(v_list[i])
        elif search_type == 1: #eSummary query term
            url_query = eSummary_generate_query(v_list[i])
            #If the variant has no IDs, skip current variant
            if url_query == 1:
                continue

        #If 3 requests have been made, check to see elasped time
        if (request_counter >= 3):
            current_time = time.time() #get current time
            time_diff = current_time - start_time #get the difference in time
            if (time_diff < 1.00): #Check to see if 1 second has elasped
                time.sleep(1.0 - time_diff) #If not, slepp the rest time
            #Reset the start time and counters
            request_counter = 0
            start_time = time.time()

        #Perform the appropriate search and process results
        if search_type==0: #Access eSearch
            #Search and retrieve the results in a List
            result = eSearch_getIDs(url_query)
            #Call function to initiate the list into the variant
            eSearch_processResults(v_list[i], result)
        elif search_type==1: #Access eSummary
            #Search and retrieve results in a dictionary
            result_dict = eSummary_getResult(url_query)
            #Check for exception and potential program termination
            if result_dict == 1:
                return 1
            #Input nested dictionary into the Var object instance
            v_list[i].recordLib = result_dict

        #Increase the request counter by 1
        request_counter += 1
    #After loop, return success
    return 0

""" #### Helper Functions to find presence of and retrieve IDs ####"""
#Function to generate the eSearch-appropriate url query from a Var object
def eSearch_generate_query(var):
    #Combine all searchable items into a single list
    """
    search_term_list = var.anno_list #Direct copy for setting the base list
    if any(char.isdigit() for char in var.snp): #if it has an rs number
        search_term_list.append(var.snp) #append the rs number
    """
    search_term_list=[var.chromosome+"[chr]", var.position+"[chrpos37]"]

    #The base url to access the ClinVar database via EUtils
    url_base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term="

    #URL-encode then join the search terms into an "OR"-separated string
    encoded_terms = "("+urllib.quote_plus(") AND (".join([t for t in search_term_list]))+")" #NOTE: AND or OR

    #Add the max return number (temp one for now)
    retmax = "&retmax=500"

    #Combine and return all the search terms
    return url_base+encoded_terms+retmax


#Function that access eSearch and return a list of ClinVar IDs
def eSearch_getIDs(url_query):
    #Access ClinVar eSearch via url and retrieve xml data
    html = urllib.urlopen(url_query).read()
    #Parse html XML data into xml tree object
    root = ET.fromstring(html)
    #Find and access the IdList element
    for IdList in root.iter("IdList"):
        #Return all of the ClinVar record IDs as a List
        return [Id.text for Id in IdList.iter("Id")]
    #Do further error processing?
    return 1

#Function to interpret and initiate the results from eSearch
def eSearch_processResults(var, result):
    #If a list is returned with elements inside
    if isinstance(result, list) and len(result)!=0:
        #Initiate the results in the Var objects
        var.IdList = result
    #If empty list, no results are found
    elif isinstance(result, list) and len(result)==0:
        #Initiate empty list in the Var objects
        var.IdList = result
    #If other cases, unexpected events happened
    else:
        print "[ERROR] Something went wrong"

"""#### Helper Functions to find summary of records and pathogenicity ####"""
#Function to generate the eSummary url query from a Var object
def eSummary_generate_query(var):
    #If the Var has no record, return appropriately
    if len(var.IdList)==0:
        return 1
    #Combine string of IDs from the variant (If the Var has record)
    id_string = ",".join([Id for Id in var.IdList])
    #Base URL to access eSummary
    url_base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id="
    #Return the url and IDs
    return url_base+id_string

#Function that access eSummary and return pathogenicity status(es)
def eSummary_getResult(url_query):
    #Access ClinVar eSummary via url and retrieve xml data
    html = urllib.urlopen(url_query).read()
    #Parse html XML data into xml tree object
    root = ET.fromstring(html)
    #Create a dictionary to store the information
    doc_set = {}
    #Iterate through each record document found
    for doc in root.iter('DocumentSummary'):
        #Create a sub-dictionary for each piece of information
        doc_sum = {}
        #Store the clnical significance
        try:
            doc_sum['clin_sig']=doc.find('clinical_significance').find('description').text
        except AttributeError as e: #If clin significance not found...
            #Print error message
            print "AttributeError ({0}): {1}. Initiated as 'None'".format(e.errno, e.strerror)
            #Initiate NoneType as attribute
            doc_sum['clin_sig'] = None
        #For all other exceptions, return message (to terminate program)
        except Exception as e:
            print e
            return 1

        #Store the condition(s) as an array
        try:
            doc_sum['cond']= [t.find('trait_name').text for t in doc.find('trait_set').iter('trait')]
        except AttributeError as e: #If condition(s) not found...
            #Print error message
            #print "AttributeError ({0}): {1}. Initiated as 'None'".format(e.errno, e.strerror)
                #Decided against returning error message since sometimes there
                #   can be empty conditions
            #Initiate NoneType as attribute
            doc_sum['cond'] = None
        #For all other exceptions, return message (to terminate program)
        except Exception as e:
            print e
            return 1

        #Input the clinical significance and conditions into document set
        doc_set[doc.get('uid')] = doc_sum
    #Return the nested dictionary
    return doc_set
