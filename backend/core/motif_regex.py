import re

COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVNacgtryswkmbdhvn', 'TGCAYRWSMKVHDBNtgcayrwsmkvhdbn')


class MotifRegex:
    def __init__(self, regex_pattern: str):
        self.regex_pattern = regex_pattern.upper()
        self.forward_regex_pattern = self._get_forward_pattern(self.regex_pattern)
        self.reverse_regex_pattern = self._get_reverse_pattern(self.regex_pattern)

    def _get_forward_pattern(self, regex_pattern: str) -> re.Pattern[str]:
        parsed_pattern = regex_pattern.replace('N', '[ACGT]')
        return re.compile(f"(?=({parsed_pattern}))", re.IGNORECASE)

    def _get_reverse_pattern(self, regex_pattern: str) -> re.Pattern[str]:
        parsed_complement_pattern = regex_pattern.replace('N', '[ACGT]').translate(COMPLEMENT)

        tokens = re.findall(r'\[.*?\]|.', parsed_complement_pattern)

        reverse_complement_pattern = "".join(reversed(tokens))

        return re.compile(f"(?=({reverse_complement_pattern}))", re.IGNORECASE)

    def __repr__(self):
        return f"MotifRegex(pattern='{self.regex_pattern}')"
