# ClinVar_PathoID
Searches up a list of variants to identify clinical significance and (if present) condition(s).

# Instruction for use:
1. Download the three .py files below into a directory of your choice
  - CV_PathoID.py
  - variant.py
  - connect.py

2. **(Optional, if gene-specific filtering is desired)** Create a plaintext file named *Wanted_Genes.txt* in the directory that the python (.py) files are in
  - In each new line of the file, enter a new gene that you wish to filter for (I.e. keep variants with this gene)
  - Example format below:
    ```
    GENE1
    GENE2
    GENE3
    ```
  - *Note that* you can create this file even if you do not wish to filter for genes. There are additional user-prompts as the program is running which will confirm whether you want to filter for genes or not.

3. Use the command line to access the directory the .py files are located in. Once you are there, run the program from command line by typing the following followed by the "enter" key:
  ```
  python CV_PathoID.py path_to_your_input_variant_file
  ```

4. Follow the prompts from the program. The output file(s) will be found in the same directory as the input files.


# Note on variant annotation formatting
Below are the steps taken (in the order they appear in) to format each detailed variant annotation into what is inputted into ClinVar to search.
- Replace ":exon" with ".", if present
- If there are brackets present in the annotation, simply return what is contained in the brackets (without the next steps)
- Remove everything before "NM" and after (inclusive) ":p"
- Change the nucleotide variation formatting:
  - I.e. ":c.G123C" is changed into ":c.123G>C"
- Example:
  - Raw annotation: *exonic:nonsynonymous_SNV:SRY:NM_003140:exon1:c.A593C:p.Y198S*
  - Formatted annotation: *NM_003140.1:c.593A>C*

In accordance with the above formatting steps, the following assumption about each of the raw variant annotation is made:
- All variants are of the "NM_" type (no "NR_", or others)
- All variants have single nucleotide variations
- Bracketed variants provide the correct format for ClinVar search

Note that the above formatting does not provide a comprehensive translation into a ClinVar-appropriate variant annotation format, and is therefore subjected to change. Feel free to to change the formatting steps in the file **variant.py** within the function **format_annotation()**.

# Other notes
File functionality descriptions:
 - CV_PathoID.py: carries out the main input/output and function calls
 - variant.py: contain the object classes and related helper functions
 - connect.py: related functions to connect to and access ClinVar

Assumptions made about each row of variants of the input file:
- I (shamefully) hardcoded the columns of interest
  - Content location:
    - 10th column contains the gene symbol
    - 8th column contains the variant function type (VFT)
    - 11th column contains the detailed variant annotation
    - 12th column contains the "All SNPs" (rs number)
  - Some things to note regarding the order:
    - The inputted column index (in the source code) is -1 from the actual column number, due to the fact that Lists index start from 0 instead of 1
    - The order of input must follow the above: gene name, function type, detailed annotation, and rs number
- Detailed variant annotation formatting
  - If there are multiple detailed annotations, each one is separated via the pipe ("|") character
  - There are no commas (",") in the detailed variant annotation, else the .csv output file will likely be compromised. This case is not seen yet but one should be careful in the future.

Other notes:
- Currently, each variant is search via a separate request (which includes all of its annotations and rs number)
- The max speed is fixed to 3 request / second, to adhere to the NCBI guideline to avoid excessive requests. Actual runtime will depend on the internet and processor speed but is likely to be close to the max runtime.
- The program will access eSearch and eSummary separately, in their respective order. eSearch is used to find (if available) a list of ClinVar record IDs (or indicate that no records are found); while eSummary will use the generated ID list to find pathogenicity status for the variant.
- If multiple records are found for a variant, all pathogenicity and disease conditions will be compiled together.
  - Format: [Patho 1]|[Patho 2]|[Patho 3] (file delimiter) [Condition 1]|[Cond 2]...
  - Note that there will be no way of knowing exactly how to distinguish between these individual records other than manually searching them
- The input filename cannot have any periods (".") other than right in front of the file extension (e.g. file.output and file.csv)
  - The file has to also be in your current directory
- No exception handling against networking / socket errors. Unsure how necessary this functionality is at this point.
- Likely error: in the URL functions: eSearch_getIDs and eSummary_getResult of connect.py
  - For now, any exceptions that arises other than AttributeError will cause the entire program to terminate immediately. May wish to fix this in the future to catch more exceptions.
