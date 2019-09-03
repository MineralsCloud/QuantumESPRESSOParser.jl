# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/qeinputparser.py
const NAMELIST_BLOCK_REGEX = r"""
^ [ \t]* &(\S+) [ \t]* $\n  # match line w/ nmlst tag; save nmlst name
(
 [\S\s]*?                # match any line non-greedily
)                        # save the group of text between nmlst
^ [ \t]* / [ \t]* $\n    # match line w/ "/" as only non-whitespace char
"""mx
