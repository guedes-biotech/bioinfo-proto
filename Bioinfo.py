from typing import Union
from tkinter import filedialog
import csv

#variables that help to build the csv file
dict_keys_from_fasta_file = list()
there_are_uniprot_files = there_are_other_database_files = global_flag_for_dict_when_full_data_function = False

#Dictionary that relate the biological information processing 
biological_dictionary = {"dna-dna" : {"A":"T", "T":"A", "C":"G", "G":"C"},
                         "dna-rna" : {"A":"U", "T":"A", "C":"G", "G":"C"},
                         "rna-dna" : {"U":"A", "A":"T", "C":"G", "G":"C"},
                         "one letter code" : {
                                                "Arg" : "R", "Asn" : "N", "Asp" : "D", "Cys" : "C",
                                                "Glu" : "E", "Gln" : "Q", "Gly" : "G", "His" : "H",
                                                "Ile" : "I", "Leu" : "L", "Lys" : "K", "Met" : "M",
                                                 "Phe" : "F", "Pro" : "P", "Ser" : "S", "Thr" : "T",
                                                "Trp" : "W", "Tyr" : "Y", "Val" : "V", "Ala" : "A"
                                             },

                         "codons" : {
                                        "Phe": ["UUU", "UUC"],
                                        "Leu": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
                                        "Ile": ["AUU", "AUC", "AUA"],
                                        "Met": ["AUG"],
                                        "Val": ["GUU", "GUC", "GUA", "GUG"],
                                        "Ser": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
                                        "Pro": ["CCU", "CCC", "CCA", "CCG"],
                                        "Thr": ["ACU", "ACC", "ACA", "ACG"],
                                        "Ala": ["GCU", "GCC", "GCA", "GCG"],
                                        "Tyr": ["UAU", "UAC"],
                                        "His": ["CAU", "CAC"],
                                        "Gln": ["CAA", "CAG"],
                                        "Asn": ["AAU", "AAC"],
                                        "Lys": ["AAA", "AAG"],
                                        "Asp": ["GAU", "GAC"],
                                        "Glu": ["GAA", "GAG"],
                                        "Cys": ["UGU", "UGC"],
                                        "Trp": ["UGG"],
                                        "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
                                        "Ser": ["AGU", "AGC"],
                                        "Gly": ["GGU", "GGC", "GGA", "GGG"],
                                        "Stop": ["UAA", "UAG", "UGA"]
                                    },

                         "amino acids full name" : {
                                                        "R" : "Arginine", "N" : "Asparagine", "D" : "Aspartic Acid",
                                                        "C" : "Cysteine", "E" : "Glutamic Acid", "Q" : "Glutamine",
                                                        "G" : "Glycine", "H" : "Histidine", "I" : "Isoleucine",
                                                        "L" : "Leucine", "K" : "Lysine", "M" : "Methionine",
                                                        "F" : "Phenylalanine", "P" : "Proline", "S" : "Serine",
                                                        "T" : "Threonine", "W" : "Tryptophan", "Y" : "Tyrosine",
                                                        "V" : "Valine", "A" : "Alanine"},
                        } 
def check_sequence_type(sequence : str):
        """ It returns if the sequence is "DNA", "RNA" or "Amino acid"

        Args:
            sequence (string): biological sequence
        """
        if any(letter in sequence for letter in "R N D Q E H I L K M F P S W Y V".split()):
            return "Amino acid"
        elif "U" in sequence:
            return "RNA"
        else:
            return "DNA"


def fasta2sequence(fasta : Union[str, list], full_data = False, doc = "") -> Union[str, dict, list, None]:
    """
    This function returns biological sequences from a fasta file. You can use arguments to get other info about the sequences and/or save the results in a csv or txt file.

    Args:
        fasta (string or list): name(s) of the fasta file(s)

        full_data (boolean): If True, returns a dictionary with the sequences, sequence codes, the species of origin, the file of origin and the type of the sequence (amino acid, RNA or DNA).

        doc (string): can be either "csv" or "txt". It will save the results in a file of the format you choose 
    """
    #Process the data from the file
    def get_data():
        """auxiliary function that gets the data from a fasta file
        """

        #Info from the header
        def from_header():
            """
            This function will get the information from the header of a sequence. Output will change based of the database you get the files,  NCBI sequences have diferent headers from UniProt sequences and other databases. 
            """
            #variables that help to build the csv file
            global dict_keys_from_fasta_file, there_are_uniprot_files, there_are_other_database_files
            
            #Headers from Uniprot have OS(Origin Species) inside them      
            if "OS=" in header:
                output = {"Sequence code" : header.split("|")[1],
                          "Species of origin" : " ".join(header.split("OS=")[1].split()[:2])
                          }
                if not there_are_uniprot_files:
                    for item in output.keys():
                        dict_keys_from_fasta_file.append(item)
                    there_are_uniprot_files = True
            #Headers from other databases usually don't have a padronized way to specify the information, this part of the code was made based on NCBI headers
            else:
                output = {"Sequence code" : header.split()[0][1:], 
                          "description" : " ".join(header.split()[1:])}
                if not there_are_other_database_files:
                    for item in output.keys():
                        dict_keys_from_fasta_file.append(item)
                    there_are_other_database_files = True
            return output
        
        #Structure of the full_data output
        def dict_when_full_data():
            """Dictionary builder
            """
            global global_flag_for_dict_when_full_data_function, dict_keys_from_fasta_file
            output = {}
            output.update(from_header())
            dict1 = {"File of origin" : fasta,
                    "Sequence type" : check_sequence_type(sequence),
                    "Sequence" : sequence
                    }
            
            #It adds the keys above to the csv headears
            if not global_flag_for_dict_when_full_data_function:
                for element in dict1.keys():
                    dict_keys_from_fasta_file.append(element)
                global_flag_for_dict_when_full_data_function = True

            output.update(dict1) 
            return output

        with open(fasta, "r") as p: #Count how many sequences have in the file
            lines = p.readlines()
            n_sequences = sum(">" in line for line in lines)
            p.seek(0)
            for line in lines:
                if line == "":
                    lines.remove(line)   

        #File has one biological sequence and I only want the sequence
        if not full_data and n_sequences == 1:
            with open(fasta, "r") as p:
                p.readline()
                output = p.read().replace("\n", "")
                
            
        #File has more than one biological sequence and I only want the sequence
        elif not full_data and n_sequences > 1:
            headers_index = []
            output = []
            for k, element in enumerate(lines): #Get the line number of all headers 
                if ">" in element:
                    headers_index.append(k)
            with open(fasta, "r") as p:
                for n in range(n_sequences): #Make one string with all the lines between two headers
                    sequence = ""
                    if n == n_sequences - 1: #The last header
                        p.readline()
                        output.append(p.read().replace("\n", ""))
                    else:
                        p.readline()
                        for k in range(0, headers_index[n+1] - headers_index[n] - 1):
                            sequence = sequence + p.readline().replace("\n", "")
                        output.append(sequence)
            

        #File has one biological sequence and I want all the data
        elif full_data and n_sequences == 1:
            with open(fasta, "r") as p:
                header = p.readline()
                sequence = p.read().replace("\n", "")
                output = dict_when_full_data()
            

        #File has more than one biological sequence and I want all the data
        elif full_data and n_sequences > 1:
            headers_index = []
            output = []
            for k, element in enumerate(lines): #Get the line number of all headers
                if ">" in element:
                    headers_index.append(k)
            with open(fasta, "r") as p:
                for n in range(n_sequences): #Make one string with all the lines between two headers
                    sequence = ""
                    
                    if n == n_sequences - 1: #Last header
                        header = p.readline()
                        sequence = p.read().replace("\n", "")
                        output.append(dict_when_full_data())
                    
                    else:
                        header = p.readline()
                        for k in range(0, headers_index[n+1] - headers_index[n] - 1):
                            sequence = sequence + p.readline().replace("\n", "")
                        output.append(dict_when_full_data())
            
        return output
    #Working with one fasta file    
    if isinstance(fasta, str):
        output = get_data()
    
    #Working with two or more fasta files
    elif isinstance(fasta, list):
        output = []
        for n in fasta:
            in_file = fasta2sequence(n, full_data= full_data)
            if isinstance(in_file, list):
                for element in in_file:
                    output.append(element)
            else:
                output.append(in_file)
    
    #In case you want to save the output in a document
    doc_type = doc.lower().strip()
    if doc_type == "txt": #txt documents
        path = filedialog.asksaveasfilename(defaultextension=".txt")
        print(f"The results below will be saved in {path}")
        with open(path, 'w') as file:
            #string output
            if isinstance(output, str):
                file.write("Sequence No.1: \n")
                file.write(f"\t{output}")
            #Dictionary output
            elif isinstance(output, dict):
                file.write(f"Sequence No.1:\n" )
                for k, s in output.items():
                    file.write(f"\t{k}: {s}\n")
            
            elif isinstance(output, list):
                #list of strings output
                if isinstance(output[0], str):
                    for k, s in enumerate(output):
                        file.write(f"Sequence No.{k + 1}:\n")
                        file.write(f"\t{s}")
                        file.write("\n"*2)
                #List of dictionarys output
                elif isinstance(output[0], dict):
                    for n in range(len(output)):
                        file.write(f"Sequence No.{n + 1}:\n")
                        for k, s in output[n].items():
                            file.write(f"\t{k}: {s}\n")
                        file.write("\n")
            
    elif doc_type == "csv": #csv documents
        path = filedialog.asksaveasfilename(defaultextension=".csv")
        print(f"The results below will be saved in {path}")
        with open(path, "w", newline="") as file:
            #String output
            if isinstance(output, str):
                writer = csv.writer(file)
                writer.writerow(["ID", "Sequence"])
                writer.writerow(["1", output])
            elif isinstance(output, dict) or isinstance(output[0], dict):
                csv_headers = []
                for element in dict_keys_from_fasta_file:
                    if element not in csv_headers:
                        csv_headers.append(element)
                csv_headers.remove("Sequence")
                csv_headers.append("Sequence")
                writer = csv.DictWriter(file, fieldnames = csv_headers)
                writer.writeheader()
                if isinstance(output, list):
                    writer.writerows(output)
                else:
                    writer.writerow(output)
            elif isinstance(output[0], str):
                data = []
                n = 1
                for element in output:
                    data.append({"ID" : n, "Sequence" : element})
                    n += 1
                writer = csv.DictWriter(file, fieldnames = ["ID", "Sequence"])
                writer.writeheader()
                writer.writerows(data)

    return output

def complementar(sequence:str) -> str:
    """This function returns the complementary DNA strand of the input sequence. It will return an empty string if it doesn't recognize the sequence as DNA or RNA.

    Args:
        sequence (str): It needs to be a DNA or a RNA sequence. When the input sequence is RNA the will function will work as the reverse transcriptase
    """
    seq_type = check_sequence_type(sequence)
    output = ""
    #It turns the output in the complementary DNA strand from the DNA sequence inputted
    if seq_type == "DNA":
        for base in sequence.upper():
            output += biological_dictionary["dna-dna"][base]
        return output
    
    #It turns the output in the complementary DNA strand from the RNA sequence inputted
    elif seq_type == "RNA":
        for base in sequence.upper():
            output += biological_dictionary["rna-dna"][base]
        return output
    
    #I'm using a print in case this function is used in a loop and I want to know which sequence is giving me trouble
    else:
            print(f"{sequence} is not a DNA or RNA sequence")
            return ""

def transcription(sequence:str) -> str:
    """This function returns the transcription of the input sequence. It will return an empty string if it doesn't recognize the sequence as DNA or RNA.

    Args:
        sequence (string): Biological sequence. It can be DNA or RNA, where RNA sequences will return the reverse transcription of the sequence
    """
    seq_type = check_sequence_type(sequence)
    output = ""
    if seq_type == "DNA":
        for base in sequence.upper():
            output += biological_dictionary["dna-rna"][base]
        return output
    elif seq_type == "RNA":
        for base in sequence.upper():
            output += biological_dictionary["rna-dna"][base]
        return output
    else:
        print(f"{sequence} is not a DNA or a RNA sequence")
        return ""
    
def translation(RNA:str, threeLcode = False) -> str:
    """This function translate a RNA sequence to an amino acid sequence

    Args:
        rna (string): RNA sequence

        threeLcode (bool): Defines the layout of the output. If it keeps False, the output will be displayed with one-letter code amino acids. otherwise, the output will be displayed with three-letter code and "-" betwen them
    """
    seq_type = check_sequence_type(RNA)
    if seq_type != "RNA":
        print(f"{RNA} is not a RNA sequence")
        return ""
    output = ""
    for n in range(0, len(RNA), 3):
        codon = RNA[n:n+3]
        if len(codon) == 3:
            for key, value in biological_dictionary["codons"].items(): 
                if codon in value:
                    output += biological_dictionary["one letter code"][key]
                    break

    if threeLcode:
        oneLcode = output
        output = ""
        for aa in oneLcode:
            for key, value in biological_dictionary["one letter code"].items():
                if aa in value:
                    output += f"-{key}"
                    break
        output = output[1:]
    return output