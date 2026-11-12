# Manuscript drafts

The Word file is the master copy of the text from V1 on (the Claude Doc draft is superseded).

Revision loop:

1. Blake reads `draft_vN.docx` in Word, adds comments (and edits, ideally with Track Changes on), and saves the file in place.
2. Claude reads the comments and their anchored text, makes the changes as tracked changes in a new file
   `draft_vN+1.docx`, and answers each comment in its thread.
3. Figures come from `scripts/figures/` (`output/figures/*.png`); when a figure changes, its image is replaced in the next draft.

Line numbers are continuous, so comments can also say "line 123".

Clean copies for co-authors (`draft_vN_clean.docx`; the `.pdf` beside it is gitignored) have every change accepted and no comments. They also
leave out the "About this draft" notes, the open-items list and inline work notes. The V5 clean copy also fills in the Methods facts (animals, surgery, display, fixation, mapping, recording), taken
from Mitchell et al. 2023 and checked against the recordings' MonkeyLogic configuration and the task code. It is the
base for V6.

V6 starts from Blake's comments on the V5 clean copy (`draft_v5_clean.docx` as he commented it), with the changes tracked on that copy.
