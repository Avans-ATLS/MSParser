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
        inchi: str = None,
        formula: str = None,
        retention_time: float = None,
        ccs: float = None,
        ionmode: str = None,
        compound_class: str = None,
        comment: str = None,
        n_peaks: int = None,
        peaks: List[Tuple[float, int]] = [],
    ):
        self.name = name
        self.precursor_mz = precursor_mz
        self.precursor_type = precursor_type
        self.smiles = smiles
        self.inchi = inchi
        self.formula = formula
        self.retention_time = retention_time
        self.ccs = ccs
        self.ionmode = ionmode
        self.compound_class = compound_class
        self.comment = comment
        self.n_peaks = n_peaks
        self.peaks = peaks if peaks else []

    def __str__(self):
        return f"Record:\n{self.name}\n{self.precursor_mz}\n{self.precursor_type}\n{self.smiles}\n{self.inchi}\n{self.formula}\n{self.retention_time}\n{self.ccs}\n{self.ionmode}\n{self.compound_class}\n{self.comment}\n{self.n_peaks}\n{self.peaks}"


class MSPdb:
    """MSP database class"""

    def __init__(self):
        self.records: List[Record] = []

    def load_file(self, filename: str) -> None:
        """Read a MSP database file and add records to class

        Args:
            filename (str): path to MSP database file
        """
        f = open(filename, "r")
        items = f.read().strip().split("\n\n")
        for item in items:  # iterate over items in database
            rec = {}
            peaks = []
            for field in item.split("\n"):  # iterate over fields in item
                if ": " in field:  # fields are separated by ': '
                    key, value = field.split(": ")
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
                    name=rec.get("NAME"),
                    precursor_mz=float(rec.get("PRECURSORMZ")),
                    precursor_type=rec.get("PRECURSORTYPE"),
                    smiles=rec.get("SMILES"),
                    inchi=rec.get("INCHIKEY"),
                    formula=rec.get("FORMULA"),
                    retention_time=float(rec.get("RETENTIONTIME")),
                    ccs=float(rec.get("CCS")) if not rec.get("CCS") == "NA" else None,
                    ionmode=rec.get("IONMODE"),
                    compound_class=rec.get("COMPOUNDCLASS"),
                    comment=rec.get("Comment"),
                    n_peaks=int(rec.get("Num Peaks")),
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
        Number of unique compound classes: {len(set([record.compound_class for record in self.records]))}
        {', '.join(set([record.compound_class for record in self.records]))}
        --------------------------------
        Number of positive modes: {len([record for record in self.records if record.ionmode == 'Positive'])}
        Number of negative modes: {len([record for record in self.records if record.ionmode == 'Negative'])}
        --------------------------------
        Number of precursor types: {len(set([record.precursor_type for record in self.records]))}
        {', '.join(set([record.precursor_type for record in self.records]))}
        --------------------------------"""
        )

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

    def write_database(self, path: str) -> None:
        with open(path, "w") as outfile:
            for record in self.records:
                outfile.write(f"NAME: {record.name}" + "\n")
                outfile.write(f"PRECURSORMZ: {record.precursor_mz}" + "\n")
                outfile.write(f"PRECURSORTYPE: {record.precursor_type}" + "\n")
                outfile.write(f"SMILES: {record.smiles}" + "\n")
                outfile.write(f"INCHIKEY: {record.inchi}" + "\n")
                outfile.write(f"FORMULA: {record.formula}" + "\n")
                outfile.write(f"RETENTIONTIME: {record.retention_time}" + "\n")
                if record.ccs:
                    outfile.write(f"CCS: {record.ccs}" + "\n")
                outfile.write(f"IONMODE: {record.ionmode}" + "\n")
                outfile.write(f"COMPOUNDCLASS: {record.compound_class}" + "\n")
                outfile.write(f"Comment: {record.comment}" + "\n")
                outfile.write(f"Num Peaks: {record.n_peaks}" + "\n")
                for peak in record.peaks:
                    outfile.write(f"{peak[0]}\t{peak[1]}\n")
                outfile.write("\n")
