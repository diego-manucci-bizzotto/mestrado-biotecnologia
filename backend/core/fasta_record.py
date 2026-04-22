import re


class FastaRecord:
    def __init__(self, header: str, sequence: str):
        self.full_header = header
        self.sequence = sequence

        id_match = re.search(r'\((.*?)\)', header)
        self.gene_id = id_match.group(1) if id_match else "Unknown"

        self.region_type = header.split(' ')[0].replace('>', '')

        coord_match = re.search(r'\[(\d+-\d+)\]', header)
        self.coords = coord_match.group(1) if coord_match else None

        self.is_reverse = "[reverse complement]" in header.lower()

        self.description = header

    def __repr__(self):
        return f"FastaRecord(id={self.gene_id}, reverse={self.is_reverse}, len={len(self.sequence)})"
