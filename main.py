import os
import argparse
from sys import argv
import pandas as pd
from rdkit import Chem

# import necessary classes
from msparser.MSPdb import MSPdb, Record
from msparser.parsers import Formula

########## possible future functions to be implemented in MSParser ##########
# blanco vs sample comparison of msp files
# filter msp database by using a file with search terms and filtered database name
# get the unique ontologies/precurzorMZ/etc.  with their counts 
# filter the MS2 spectra by only keeping the peaks haveng at least one percent of the highest intensity in the spectrum


### add custom functions below here ###

def argumement_parser():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="MSParser: A tool for parsing and filtering msp files.")
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the msp file to be processed.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file name with full path.')
    parser.add_argument('-F', '--filter', type=str, required=False, help='Element to filter on (options; name, ontology, precurzor_mz, ionmode, retention_time, formula, smiles, inchikey, comment, collision_energy, number of peaks or peak intensity).')
    parser.add_argument('-t', '--threshold', type=float, required=False, help='Threshold value required if the filtering is selected.')
    parser.add_argument('-s', '--searchterm', type=str, required=False, help='searchterm used in te filtering.')
    parser.add_argument('-r', '--remove', type=str, required=False, help='Element to remove from the database (options; name, ontology, precurzor_mz, ionmode, retention_time, formula, smiles, inchikey, comment, collision_energy, number_of_peaks or peak_intensity).')
    
    return parser.parse_args()


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

#FEMKE_K_BLANCO_VS_SAMPLE
def process_files(directory, blanco_pattern='_BL_'):
    """Process files and create dictionaries with error handling.
    function to use for a directory with both blanco and sample files"""
    blanco_dict = {}
    sample_dict = {}
    
    for filename in os.listdir(directory):
        if not filename.endswith('.msp'):
            continue
            
        file_path = os.path.join(directory, filename)
        try:
            db = MSPdb()
            db.load_file(file_path)
            
            if blanco_pattern in filename:
                blanco_dict[filename] = db
            else:
                sample_dict[filename] = db
                
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
    
    return blanco_dict, sample_dict


#FEMKE_K_BLANCO_VS_SAMPLE
def count_occurrences(matches_dict):
    """Count occurrences of sample compounds."""
    counts = {}
    for blanco_compound, sample_compounds in matches_dict.items():
        for sample_compound in sample_compounds:
            counts[sample_compound] = counts.get(sample_compound, 0) + 1
    return counts

#FEMKE_K_BLANCO_VS_SAMPLE
def compare_databases(blanco_db, sample_db, rt_tolerance=0.5, mz_decimal_places=2):
    """Compare databases and track hits per sample."""
    matches = {}
    blanco_occurrences = {}
    
    for blanco_compound in blanco_db.records:
        if blanco_compound.n_peaks < 5:
            continue
            
        blanco_mz = round(float(blanco_compound.precursor_mz), mz_decimal_places)
        blanco_rt = float(blanco_compound.retention_time)
        
        matching_samples = []
        for sample_compound in sample_db.records:
            sample_mz = round(float(sample_compound.precursor_mz), mz_decimal_places)
            sample_rt = float(sample_compound.retention_time)
            
            if (blanco_mz == sample_mz and 
                abs(blanco_rt - sample_rt) <= rt_tolerance):
                matching_samples.append(sample_compound.name)
        
        if matching_samples:
            matches[blanco_compound.name] = matching_samples
            blanco_occurrences[blanco_compound.name] = len(matching_samples)
    
    return matches, blanco_occurrences


#FEMKE_K_BLANCO_VS_SAMPLE
def write_results(blanco_dict, sample_dict, output_file):
    """Write results to file with total hit tracking."""
    with open(output_file, 'w') as f:
        f.write("Blanco compound hits across all samples\n")
        f.write("--------------------------------------------------\n")
        
        total_occurrences = 0
        blanco_occurrences = {}
        sample_files = list(sample_dict.keys())
        
        for blanco_filename, blanco_db in blanco_dict.items():
            for sample_filename, sample_db in sample_dict.items():
                matches, blanco_hits = compare_databases(blanco_db, sample_db)
                blanco_occurrences.update(blanco_hits)
                total_occurrences += sum(blanco_hits.values())
                
                for blanco_compound, sample_compounds in matches.items():
                    # Format line with total hits and sample list
                    formatted_line = (
                        f'{blanco_compound} has {len(sample_compounds)} hit(s) '
                        f'in {"|".join(sample_files)}'
                    )
                    f.write(f"{formatted_line}\n")
        
        # Add summary statistics
        f.write("\nSummary Statistics:\n")
        f.write("------------------\n")
        for blanco, count in sorted(blanco_occurrences.items(), 
                                  key=lambda x: x[1], 
                                  reverse=True):
            f.write(f"{blanco}: {count} hits across all samples\n")
        
        print(f"Total occurrences across all samples: {total_occurrences}")
        print("\nSummary Statistics:")
        for blanco, count in sorted(blanco_occurrences.items(), 
                                  key=lambda x: x[1], 
                                  reverse=True):
            print(f"{blanco}: {count} hits")


###### MAIN FUNCTION CODE AREA ######       

args = argumement_parser()
# load the database
db = MSPdb()
db.load_file(args.file)

print(f"The original database summary:")
# db.summary()

filtered_db = MSPdb()
filtered_db.records = [x for x in db.clean_peaks_on_intensity(threshold=args.threshold)]
print(f"The database summary after filtering on peak intensity: \n")
filtered_db.summary()

filtered_db.write_database(args.output)



# file_location = '/home/daan/databases/stefano/List.txt'
# list_of_searchterms = parse_search_words(file_location)
##################################################################################
#FEMKE_K_BLANCO_VS_SAMPLE
#Compare femke's blanco databases with their sample databases
# blanco_files, sample_files = process_files('/home/daan/databases/femke_k/blancoVSlocation/')

# for sample_filename, sample_db in sample_files.items():
#     # create a filtered_db object for storage of the filtered records
#     filtered_db = MSPdb()
#     # loop over every compound in the sample database
#     for sample_compound in sample_db.records:
#         if sample_compound.n_peaks < 5:
#             continue
#         else:
#             # rename the compound to its mz and retention time
#             sample_compound.name = f"MZ={sample_compound.precursor_mz}|RT={sample_compound.retention_time}"
#             # empty the comment line
#             sample_compound.comment = ""
#             # append the sample compound to the filtered database
#             filtered_db.records.append(sample_compound)
#     # write the filtered database to a file
#     filtered_db.write_database(f"/home/daan/databases/femke_k/blancoVSlocation/name_comment_cleaned_files/2025APR24_{sample_filename}_NameCommentClean.msp")


# write_results(blanco_files, sample_files, 
#              '/home/daan/databases/femke_k/blancoVSlocation/2025APR23_blancoVsample_compounds.txt')

###################################################################################





###################################################################################

# for sample_name, counts in sample_occurence_counts.items():
#     print(f"{sample_name}: {counts} occurences")

# blanco_comparisons = blancohits_in_samples(blanco_db, sample_db, rt_tolerance=0.5, mz_decimal_places=2)

# for blanco_compound, sample_compounds in blanco_comparisons.items():
#     print(f"Blanco compound: {blanco_compound}")
#     print(f"Sample compounds: {', '.join(sample_compounds)}")
#     print("")


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