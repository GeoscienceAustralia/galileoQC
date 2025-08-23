(linenumbers-target)=
# Flight Line Numbering

Airborne geophysical surveys are almost always flown as a series of line segments and each segment is uniquely numbered. The numbering is used in planning, acquisition, processing and QC. The method of numbering the line varies across survey providers but is always designed to encode useful information about each segment.

Surveys flying follows planned lines and the planned lines are also numbered. The segment number always includes the planned line number.

Occasionally, a segment, or part of a segment will need to be re-flown because of problems with the data or the aircraft track. This may need doing more than once. The segment number includes a re-flight number to indicate which re-flight it is.

Most survey consist of traverse lines flown parallel to each other, and control lines (sometimes called tie-lines) flown perpendicular (or occasionally at some angle close to a right angle) to the traverse lines. The control lines are used to level the data and sometimes in QC. They are typically much more widely spaced than the traverse lines. Usually a part of the planned line number, and the corresponding part of the segment number will indicate whether the line is a traverse line or a control line.

Some surveys are divided into blocks. This division might be because of spatial separation, because parts of the total area are to be flown in different traverse directions, because part of the survey are to be flown by different aircraft, because lines were added to the survey plan after acquisition had commenced, or for some other reason. The planned and segment numbers might include a numerical code to indicate which block the line is from.

Sometimes a survey area has a hole in the middle. It might be most efficient flying to simply fly the aircraft across the gap but simply to delete the data over the gap. Logistics might mean that the planned line number is the same on both sides of the gap. A segment number might be used to indicate which side of the gap.

Doubtless there are other pieces of information that might be encoded into the line numbering.

The airborne geophysics industry has no common approach to encoding information in line numbers so generally, every survey provider will use their own system. Sometimes the method will be different for different surveys flown by the same provider. Even during a single survey, circumstances might lead to a change in line numbering system. Some geophysical software requires that line numbers be integers, and some allows decimal line numbers and this difference leads to differences in survey providers approach to line numbering.

For QC, the most critical requirement is to be able to identify the planned line, given the flown line segment. It is also useful for some QC functions to be able to identify a segment as a traverse or control line.

The `updateLineAttributes` function recognises several encoding formats for line numbering:

**ARK**

For the Bonaparte ArkEx FTG survey. Details TBD.

**SGL_NSW**
: Line Number convention "LLLL.SR" :

- 'LLLL' is the planned line number (4 digits: traverse lines, 3 digits: control/tie lines). Test lines have numbers
  greater than 6999.
- "S" represents the pre-planned segment number of the line. A planned line that is not continuous will have
      separate segments with different segment numbers.
- "R" represents the particular occasion that a line/line segment was flown, starting from zero.
      For example, 1245.20 would be the planned segment '2' of the planned line 1245 flown for the first time.
      If this line segment were flown as two partial lines, their identifying line numbers would be 1245.20 and 1245.21.

**SGL_GA**

: Line Number convention "LLLL.SR" :

- 'LLLL' is the planned line number (4 digits: traverse lines, 3 digits: control/tie lines). Test lines have numbers
  greater than 6999.
- "S" represents the preplanned segment number of the line. A planned line that is not continuous will have separate segments with different segment numbers. Segment numbers are used to identify survey blocks.
- "R" represents the particular occasion that a line/line segment was flown, starting from zero.

**Xcal_nsw**
: Line Number convention "LLLLSR" :

- If the second digit is "9", then the line is a control, otherwise a traverse.
- The planned line number is "LLLL0".
- The last digit "R" is the re-flight number.
- The second last number "S" is the segment number.
- If the line number has only 3 digits, then it is a test line.

**Xcal_can**
: Line Number convention "LLLLR" :

- If the second digit is "9", then the line is a control, otherwise a traverse.
- The planned line number is "LLLL0".
- The last digit "R" is the re-flight number.

**NRG**
: Line Number convention "LLLLR" :

- The first digit is "1" for traverses, "9" for control lines.
- The last digit is the re-flight number.
- The planned line number is the same as the flown line number with re-flight number "0".



