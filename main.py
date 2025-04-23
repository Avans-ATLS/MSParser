import os
import argparse
from sys import argv
import pandas as pd

# from rdkit import Chem

from msparser.MSPdb import MSPdb, Record
from msparser.parsers import Formula

### add custom functions below here ###

def add_13c_headlabel(records: list[Record]) -> list[Record]:
    """Add a choline group to all records in a list

    Args:
        records (List[Record]): list of records
    """
    new_records = []
    for record in records:
        record.name = f'hl{record.name}'
        record.precursor_mz = record.precursor_mz + 2.00671
        record.smiles = record.smiles.replace("CC[N+]", "[13C][13C][N+]")
        molecule = Chem.MolFromSmiles(record.smiles)
        record.inchikey = Chem.MolToInchiKey(molecule)
        formula = Formula(record.formula)
        formula.change_element_count('C', '-', 2)
        formula.change_element_count('C', '+', 2, isotope=13)
        record.formula = formula.to_string()
        record.peaks = [(p[0] + 2.00671, p[1]) for p in record.peaks]
        new_records.append(record)
    
    return new_records


def unique_ontologies(records: list[Record]):
    """return a list of all unique ontologies found in a .msp database
    
    records: list of record objects from the .msp database file
    """
    individual_ontologies = []
    for record in records:
        if not record.ontology in individual_ontologies:
            individual_ontologies.append(record.ontology)

    return individual_ontologies

def update_comment_samplelocations(peak_sample_file, records: list[Record]):
    """update the comment section with sample locations

    args:
        peak_sample: txt file with peak id and the sample locations it was found in
        records: list of record objects from the .msp database file
    """
    try:  
        with open(peak_sample_file, 'r') as sample_content:
            lines = sample_content.readlines()
        
        # dictionary to store all sample iDs
        peak_samples = {}

        # process each line seperately
        for line in lines:
            
            #skip empty lines
            if not line.strip():
                continue

            # validate that the line has at least one value
            line_parts = line.strip().split('\t')
            # print(line_parts)
            if not line_parts:
                continue

            compound_ID = line_parts[0]

            if compound_ID not in peak_samples:
                peak_samples[compound_ID] = set()
            peak_samples[compound_ID].update(line_parts[1:])


        # Process each record and update its comment
        for record in records:
            try:
                # Extract current compound ID from comment
                comment = record.comment.strip()
                
                # Validate comment format
                comment_parts = comment.split('|')
                if len(comment_parts) < 2:
                    print(f"Warning: Invalid comment format: {comment}")
                    continue
                
                id_part = comment_parts[1].split('=')
                if len(id_part) != 2 or id_part[0] != 'PEAKID':
                    print(f"Warning: Invalid ID format: {comment}")
                    continue
                
                compound_id = id_part[1]
                
                # Update comment with all sample locations for this peak ID
                if compound_id in peak_samples:
                    new_comment = f"PEAKID={compound_id}|{'|'.join(sorted(peak_samples[compound_id]))}"
                    record.comment = new_comment
                    print(f"Updated record: {record.comment}")
            
            except (IndexError, AttributeError) as e:
                    print(f"Error processing record: {e}")
                    continue
        
        print("Peak ID updates completed successfully")
        return records
    
    except FileNotFoundError:
        print(f"Error: Could not find file {peak_sample_file}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")   

def parse_search_words(file_path):
    """extract the words needed for your filtering steps as a list

    args:
        file_path: full path to the file 
    """       
    with open(file_path, 'r') as of:
        terms = []
        for line in of.readlines(): # turn every line into a list and loop through the list
            line = line.strip('", \n')
            terms.append(line)

    return terms

def blanco_vs_location(blanco_db, sample_db):
    """compare the blanco database with the sample database

    args:
        blanco_db: MSPdb object with the blanco database
        sample_db: MSPdb object with the sample database
    """
    blanco_vs_sample_dict = {}
    for compound in blanco_db.records:
        # loop over the sample database and compare precurzor mz (rounded to 2 decimals) and Retention Time (rounded to 2 decimals)
        for sample_compound in sample_db.records:
            

###### MAIN FUNCTION CODE AREA ######       
# open the reference database and create the database object


# file_location = '/home/daan/databases/stefano/List.txt'
# list_of_searchterms = parse_search_words(file_location)
##################################################################################
# Check if the blanco compounds are present in the sample database, based on the PRECURSORMZ and RT
blanco_db = MSPdb()
blanco_db.load_file('/home/daan/databases/femke_k/blancoVSlocation/Msp_BL_MQ_1.msp')

sample_db = MSPdb()
sample_db.load_file('/home/daan/databases/femke_k/blancoVSlocation/Msp_G01_09Jan25.msp')




##################################################################################
# create a list of the unique ontologies in the database
# ontologies = unique_ontologies(db.records)

# db_femke = MSPdb()
# db_femke.load_file('/home/daan/databases/femke_k/Msp_2025_03_24_08_50_53_AlignmentResult_2025_03_11_10_01_19.msp')
# # db_femke.summary()

# updated_db_femke = MSPdb()

# sample_locations_percompound = '/home/daan/databases/femke_k/peak_id_per_samplecontent.txt'
# updated_records = update_comment_samplelocations(peak_sample_file=sample_locations_percompound, records=db_femke.records)

# for record in updated_records:
#     updated_db_femke.records.append(record)

# updated_db_femke.write_database('/home/daan/databases/femke_k/location_linked_compounds.msp')
##################################################################################



###### filter the ref database with your criteria ######

#TODO get the filters not hard coded

# #create a new database object for your filters
# filtered = MSPdb()

# for ontology_term in list_of_searchterms:
#     # append the records with these search terms into the filtered db
#     for record in db.filter_ontology(ontology_term):
#         filtered.records.append(record)

# for record in db.filter_ontology("Acyloins"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Beta carbolines"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Beta hydroxy acids and derivatives"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Cembranolides"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Hormada alkaloids"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Hydroxy fatty acids"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Jasmonic acid"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Methionine and derivatives"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Methyl-branched fatty acids"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Organic acids"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Oxidized fatty acids"):
#     filtered.records.append(record)
# for record in db.filter_ontology("xylenes"):
#     filtered.records.append(record)
# for record in db.filter_ontology("Sulfuric acid monoesters"):
#     filtered.records.append(record)



# filtered.summary()

# db.records.extend(filtered)
# db.summary()

# filtered.write_database('/home/daan/databases/stefano/glycosides.msp')	