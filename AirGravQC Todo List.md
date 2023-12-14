# AirGravQC Todo List

0. Get everything setup as a local package under pip.
2. Change all QC routines to report in the same format (Projectname, Blockname, summary [, detail text [, plot]]).
3. Standardise plots as much as possible (plot titles, axis labels, use of units). Write and use common plot functions.
4. Design geoWhizz v1.1 (based on `gspy`), and then upgrade to it.
5. Test the various height checks: checkSafeClearance, checkClearance, checkVert; that they work as intended.
7. `checkGNSS` only finds at most one error per line. It should find and report against the actual spec like `checkVertPlan`. Similarly it should have an option to plot lines with failures; and an option to test a specified array of lines.
8. All `check***` routines to have options `lines=[]`, `plot_flag=False`, `known=''`, `verbose=False`.
9. Standardise QC (check) functions to report: 1 - what it is doing; 2 - a one sentence result summary; 3 - line-by-line summary.
10. Get permissions from Andy Gabell to use some field data for tutorial.
17. Write checkLineSeparation() against separation = nominal +/- allowance for distance > max_distance.
18. plotBoxWhisker - x axis formatter to include decimals, write another formatter function to do this.
19. PROGRESS. All the source files are way too big - re-factor!!
20. checkXYPlan - gets confused about number of zeros after decimal point in line number
22. If the XYZ file contains a date in the format YYYY/MM/DD, then convert it to a decimal date string. See the stub `translate_date` in `gw`.
23. Write a function to compare two whizz datafiles and report which lines have differences. Include a `detail` flag, when true, print the first exemplar difference on the line.
24. checkGaps() - modify to allow gaps smaller than some minimum size.
25. commonErsHdrErrors() - include a check that a Units field is present.
26. Set up as a package on Github.
27. checkHeading plots output from one line even when there are 0 failures.
28. Have `asegToHDF` read string fields such as dates '12/04/2021' and store them.
29. Have `asegToHDF` use the key attributes to populate CoordFrame.
30. `checkSpikes` needs to give some report when there are no spikes found.
31. `checkDiurnal` needs a report when successful.

1. DONE. Get the ASEG-GDF to Whizz converter working well. (Currently reads data but ought to automatically import meta-data as well.)
6. DONE. Is `checkDrape` useful or redundant? Required for the case where the planned `drape` is included as a channel in the db.
11. DONE. Stop using "erm" - use "grd" instead.
12. DONE - Get permissions from Xcalibur to use Canobie field data for tutorial.
13. DONE - Get permissions from SGL to use GA/GSV field data for tutorial.
14. DONE - Get permissions from Bell Geospace to use some field data for tutorial.
15. PROGRESS. Make code insensitive to the case of the channel names.
16. DONE - Make flightChannel default = attrib.flight
21. DONE. Everywhere: replace 'pf' and 'mhd' as shortcuts for 'pointfiles' with 'gw'.



## geoWhizz v1.1

0. Use `gspy` as the basis.
1. Get rid of separate CoordinateFrame and hold that info as project metadata.
2. Put in EPSG metadata or other Coordinate Reference Systems (CRS) as per `gspy`.
3. Allow user to invent and store their own metadata at Project, Line, and Channel level. In other words, some metadata is there "just for the record" and is not used by the s/w and we can allow anything the user might want (and of course, `reportWhizz` will report it by design), and some metadata is actually used by the s/w and must follow the implicit conventions.
4. Stop storing channels that are expected, and are constant: `FLIGHT`, `LINE`, `PROJECTNUM`, `DOY`, etc. Channels constant across a line should be Line metadata and those constant across the project should be Project metadata. In every case, check that they are constant before reducing to metadata; report the error if they are not and do not convert to metadata in that case.
5. Channels that should be integers (e.g. observed number of GNSS satellites) must be checked that they are integers by `xyzToHDF`, and stored as integers.
6. Convert dates (e.g. DOY, decimal year, etc) to ISO date using Python's `DateTime` using function `doyToISO8601`.

## Other

13. If there are lots of lines, and an error persists over all (or most) lines, then the report can be very long. Provide an `verbose=False` parameter which, if true, limits the report to the summary.
14. Calculating the statistics of angular variables has problems with 0, 2pi. Find a fix and implement.
16. Automate setting of channel attributes by reading a channel description file.
17. The Vinton Dome database does not contain flight numbers. reportFLights() BUG - should check for 'FLIGHT' attribute and report clean result instead of crashing. gw.reportFlights(dh, detailed=True)
1. DONE? Vertical axis of second plot needs to be in `xxx xxx m` format.
2. DONE. Ditto for checkDrape, checkClearance.
3. DONE. Ditto for checkSpeed.
4. DONE. Ditto for checkRepeatLines.
5. DONE. checkRepeats to be re-factored with the plot function (called four times) separate.
6. DONE. checkSpeed plots need horiz axis to match criteria (maxDistance -> distance; maxDuration -> time).
7. DONE. All QC routines plotting a box-whisker plot to do it via a call to the same wp.plotBoxWhisker() function.
8. DONE. plotBoxWhisker to have a label for the x axis.
9. DONE. diffNoiseVturb to have a label_lines Bool parameter which labels each point on the plot so that one can identify outliers.
10. DONE. ilsNoiseVturb to have a label_lines Bool parameter which labels each point on the plot so that one can identify outliers.
11. DONE. Have `allChanStats` capable of reporting statistics on channels after mean removal and/or first differencing.
12. DONE. checkSafeClearance needs a plot_flag parameter.
15. DONE - `checkRawAGG` references the Lawin QC report. Replace this with an actual description in the Methods documentation.
18. DONE. `checkOverlaps` The report could be more informative here - "all overlaps met the requirement", say
