from typing import List, Tuple, Generator
from dataclasses import dataclass


class Record:
    """Record class for MSP database"""

    def __init__(
        self,
        name: str = None,
        precursor_mz: float = None,
        precursor_type: str = None,
        smiles: str = None,
        inchikey: str = None,
        formula: str = None,
        ontology: str = None,
        retention_time: float = None,
        ccs: float = None,
        ionmode: str = None,
        comment: str = None,
        n_peaks: int = None,
        peaks: List[Tuple[float, int]] = None,
    ):
        self.name = name
        self.precursor_mz = precursor_mz
        self.precursor_type = precursor_type
        self.smiles = smiles
        self.inchikey = inchikey
        self.formula = formula
        self.ontology = ontology
        self.retention_time = retention_time
        self.ccs = ccs
        self.ionmode = ionmode
        self.comment = comment
        self.n_peaks = n_peaks
        self.peaks = peaks if peaks is not None else []

    def __str__(self):
        return f"Record:\n{self.name}\n{self.precursor_mz}\n{self.precursor_type}\n{self.smiles}\n{self.inchikey}\n{self.formula}\n{self.ontology}\n{self.retention_time}\n{self.ccs}\n{self.ionmode}\n{self.compound_class}\n{self.comment}\n{self.n_peaks}\n{self.peaks}"


class MSPdb:
    """MSP database class"""

    def __init__(self):
        self.records: List[Record] = []

    def load_file(self, filename: str) -> None:
        """Read a MSP database file and add records to class

        Args:
            filename (str): path to MSP database file
        """
        with open(filename, "r") as f:
            items = f.read().strip().split("\n\n") #split records by blank line

        for item in items:  # iterate over items in database
            rec = {}
            peaks = []
            
            for field in item.split("\n"):  # iterate over fields in item
                if ": " in field:  # fields are separated by ': '
                    key, value = field.split(": ", 1)
                    key = key.strip()
                    value = value.strip()
                    rec[key] = value

                elif "\t" in field:  # peaks are separated by tab
                    peaks.append(tuple(map(float, field.split("\t"))))

                if "CCS" not in rec:  # some records do not have CCS
                    rec["CCS"] = "NA"

            # create record from rec and peaks
            try:
                record = Record(
                    name=rec.get("NAME", ""),
                    precursor_mz=float(rec.get("PRECURSORMZ", 0)),
                    precursor_type=rec.get("PRECURSORTYPE", ""),
                    smiles=rec.get("SMILES", ""),
                    inchikey=rec.get("INCHIKEY", ""),
                    formula=rec.get("FORMULA", ""),
                    ontology=rec.get("ONTOLOGY", ""),
                    retention_time=float(rec.get("RETENTIONTIME", 0)),
                    ccs=float(rec.get("CCS")) if not rec.get("CCS") == "NA" else None,
                    ionmode=rec.get("IONMODE"),
                    comment=rec.get("COMMENT", ""),
                    n_peaks=int(rec.get("Num Peaks", 0)),
                    peaks=peaks,
                )
                self.records.append(record)
                f.close()
            except:
                print("Record failed:", item.split("\n"))  # if invalid, print record
                f.close()

    def summary(self):
        """Print a summary of the MSP database"""
        print(
            f"""
        MSP Database Summary
        --------------------------------
        Number of records: {len(self.records)}
        --------------------------------
        Number of records with CCS: {len([record for record in self.records if record.ccs is not None])}
        --------------------------------
        Number of positive modes: {len([record for record in self.records if record.ionmode == 'Positive'])}
        Number of negative modes: {len([record for record in self.records if record.ionmode == 'Negative'])}
        --------------------------------
        Number of precursor types: {len(set([record.precursor_type for record in self.records]))}
        {', '.join(set([record.precursor_type for record in self.records]))}
        -------------------------------- 
        Number of unique ontologies: {len(set([record.ontology for record in self.records if record.ontology.strip()]))}
        {', '.join(set([record.ontology for record in self.records if record.ontology.strip()]))} 
        --------------------------------\n""")

    def filter_class( 
        self, compound_class: str, records: List[Record] = None
    ) -> Generator[Record, None, None]:
        """Filter compound class by exact match. Matching records are yielded.

        Args:
            compound_class (str): class to filter on.
            records (List[Record], optional): A list of records to filter. If no list is supplied, the internal database will be filtered and returned. Defaults to None.

        Yields:
            Generator[Record, None, None]: Records passing the filter.
        """
        for record in self.records if not records else records:
            if record.compound_class == compound_class:
                yield record

    def filter_ionmode(
        self, ionmode: str, records: List[Record] = None
    ) -> Generator[Record, None, None]:
        """Filter ion mode by exact match. Matching records are yielded.

        Args:
            ionmode (str): ionmode to filter
            records (List[Record], optional): A list of records to filter. If no list is supplied, the internal database will be filtered and returned. Defaults to None.

        Yields:
            Generator[Record, None, None]: Records passing the filter.
        """
        for record in self.records if not records else records:
            if record.ionmode == ionmode:
                yield record

    def unique_ontologies(self, records: list[Record]) -> dict[str, int]:
        """return a dictionary of all unique ontologies and their counts found in the database.
    
        Args:
            records: list of record objects from the .msp database file

        Returns:
            ontology_counts: dictionary of unique ontologies and their counts
        """ 
        ontology_counts = {}
        for record in self.records:  
            if not record.ontology in ontology_counts.keys():
                ontology_counts[record.ontology] = 1
            else:
                ontology_counts[record.ontology] += 1
        
        return ontology_counts  

    def filter_ontology(
        self, ontology: str, records: List[Record] = None
    ) -> Generator[Record, None, None]:
        """Filter ontology by exact match. Matching records are yielded.

        Args:
            ontology (str): ontology to filter
            records (List[Record], optional): A list of records to filter. If no list is supplied, the internal database will be filtered and returned. Defaults to None.

        Yields:
            Generator[Record, None, None]: Records passing the filter.
        """
        for record in self.records if not records else records:
            if record.ontology == ontology:
                yield record

    def ffilter_name(
        self, name: str, records: List[Record] = None
    ) -> Generator[Record, None, None]:
        """Fuzzy filter name by substring. Matching records are yielded.

        Args:
            name (str): substring to filter on.
            records (List[Record], optional): A list of records to filter. If no list is supplied, the internal database will be filtered and returned. Defaults to None.

        Yields:
            Generator[Record, None, None]: Records passing the filter.
        """
        for record in self.records if not records else records:
            if name in record.name:
                yield record

    def filter_precursor_type(
        self, precursor_type: str, records: List[Record] = None
    ) -> Generator[Record, None, None]:
        """filter records on precursor type by exact match. Matching records are yielded.

        Args:
            precursor_type (str): Precursor type to filter on.
            records (List[Record], optional): A list of records to filter. If no list is supplied, the internal database will be filtered and returned. Defaults to None.

        Yields:
            Generator[Record, None, None]: Records passing the filter.
        """
        for record in self.records if not records else records:
            if record.precursor_type == precursor_type:
                yield record

    def clean_peaks_on_intensity(
        self, records: List[Record] = None, threshold: float = 0.01
    ) -> Generator[Record, None, None]:
        """Filter peaks on intensity. Peaks below the percentage threshold are removed.

        Args:
            records (List[Record], optional): A list of records to filter. If no list is supplied, the internal database will be filtered and returned. Defaults to None.
            threshold (float, optional): Intensity threshold relative to the highest peak. Defaults to 0.01 = 1%

        Yields:
            Generator[Record, None, None]: Records passing the filter.
        """
        for record in self.records if not records else records:
            if record.peaks:
                max_intensity = max([peak[1] for peak in record.peaks])
                record.peaks = [peak for peak in record.peaks if peak[1] >= max_intensity * threshold]
                record.n_peaks = len(record.peaks)
                yield record


    def write_database(self, path: str) -> None:
        with open(path, "w") as outfile:
            for record in self.records:
                outfile.write(f"NAME: {record.name}" + "\n")
                outfile.write(f"PRECURSORMZ: {record.precursor_mz}" + "\n")
                outfile.write(f"PRECURSORTYPE: {record.precursor_type}" + "\n")
                outfile.write(f"SMILES: {record.smiles}" + "\n")
                outfile.write(f"INCHIKEY: {record.inchikey}" + "\n")
                outfile.write(f"FORMULA: {record.formula}" + "\n")
                outfile.write(f"ONTOLOGY: {record.ontology}" + "\n")
                outfile.write(f"RETENTIONTIME: {record.retention_time}" + "\n")
                if record.ccs:
                    outfile.write(f"CCS: {record.ccs}" + "\n")
                outfile.write(f"IONMODE: {record.ionmode}" + "\n")
                outfile.write(f"COMMENT: {record.comment}" + "\n")
                outfile.write(f"Num Peaks: {record.n_peaks}" + "\n")
                for peak in record.peaks:
                    outfile.write(f"{peak[0]}\t{peak[1]}\n")
                outfile.write("\n")
