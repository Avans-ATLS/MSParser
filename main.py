import os
import argparse

from rdkit import Chem

from msparser.MSPdb import MSPdb, Record
from msparser.parsers import Formula



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
        record.peaks = [(p[0] + 2.0671, p[1]) for p in record.peaks]
        new_records.append(record)
    
    return new_records
        
        
        
db = MSPdb()
db.load_file('tmpdata/MSDIAL_plusPT_pos.msp')
db.summary()

filtered = MSPdb()
for record in db.filter_class('PC'):
    filtered.records.append(record)
for record in db.filter_class('EtherPC'):
    filtered.records.append(record)
for record in db.filter_class('LPC'):
    filtered.records.append(record)
for record in db.filter_class('EtherLPC'):
    filtered.records.append(record)
for record in db.filter_class('SM'):
    filtered.records.append(record)
    
filtered.summary()

labeled_records = add_13c_headlabel(filtered.records)

db.records.extend(labeled_records)

db.summary()
db.write_database('tmpdata/MSDIAL_plusPT_pos_13C.msp')	