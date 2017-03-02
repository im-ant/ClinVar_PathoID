#!/usr/bin/python

"""##################################################################
# Helper class that contains the Var (variant) object class; also
#   carries out variant-related functions (e.g. formatting)
#
# Note to self:
#   - For now, only the detailed annotation and snp will be searched;
#       the rest will be print directly to output without any
#       intermediate manipualtions
#       - This means that the genes will not have to be split
#   - Probably should combine IdList and recordLib into a single variable
#       since recordLib is just an extension of IdList, hmmmmmmmm
#   - Not sure how to handle non-existing information workflow. Keep it as
#       a NoneType, or make user manually check?
#
#
##################################################################"""

class Var:
    #Class constructor and attributes
    def __init__(self, gene, func_type, anno, snp):
        #Content from file
        self.gene = gene #Gene symbol
        self.function_type = func_type #Function type
        self.annotation = anno #Detailed variant annotation
        self.snp = snp #All SNPs
        #Pre-search program-created contents
        self.searchable = True #Whether it is a wanted variant (yes for now)
        self.anno_list = None #Formatted list of annotation(s)
        #Post eSearch attributes
        self.IdList = None
        #Post eSummary attributes
        self.recordLib = None
    ########## Object methods ##########
    #Method to output the clinical significance
    def output_clin_sig(self):
        #If there are no clinical significance, output No items found
        if self.recordLib == None:
            return ""
        #Join all clinical significance, separated by square brackets
        #   Also replace using colons to prevent comma-separation
        clinSig = ']['.join([info['clin_sig'] for key,info in self.recordLib.iteritems()])
        return ("[%s]"%(clinSig)).replace('][',']|[').replace(',',';')

    #Method to output the condition
    def output_conditions(self):
        #If there are no conditions found, output an empty string
        if self.recordLib == None:
            return ""
        #Join all conditions together, separated by square brackets
        output_cond = ""
        for key, info in self.recordLib.iteritems():
            output_cond += "[%s]"%(']['.join(c for c in info['cond']))
        #Return the conditions, with additional pipe ('|') delimeter
        #   Also replace using colons to prevent comma-separation
        return output_cond.replace('][',']|[').replace(',',';')


#Goes through the entire variant list for pre-search formatting
def format_variantList(variant_list):
    for v in variant_list:
        #Split (potential) multiple genes by the "|" delimiter
        v.anno_list = v.annotation.split("|")
        #Iterate through each annotation to format them properly for search
        for i in range(0, len(v.anno_list)):
            v.anno_list[i] = format_annotation(v.anno_list[i])

#Function that formats an individual annotations
def format_annotation(raw_anno):
    #Temp formatting:
    #formatted_anno = raw_anno
    formatted_anno = raw_anno[raw_anno.index("NM"):]
    return formatted_anno #testline
