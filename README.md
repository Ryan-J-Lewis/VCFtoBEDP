# VCFtoBEDP
C++ code to convert the vcf.gz file format into bedp file format

vcf.gz requirements: 
1. Haplotype data must be phased (haplotype seperator = "|", not = "/") 
2. No missing haplotype data allowed (no option in bedp format)

Note: only "GT" SNP data will be converted. All others will disgaurded and printed out on the command line during execution as they are parsed. 

Instructions:
1. download the zip file of the repository and unzip in local directory
2. navigate to the .../VCFtoBEDp/VCFtoBEDp folder and enter the command "make". An exectuable, "VCFtoBEDp", will be produced and made availible in the same folder.
3. to run the program, while in the same folder, enter the command "./VCFtoBEDp", you should see the options screen print out.
4. Example of use: "./VCFtoBEDp -i input.vcf.gz -o output -O output" <br/>This command, given the input.vcf.gz file, will produce a folder named output whose content will consist of three files "output.bed", "output.bim", and "output.fam"

