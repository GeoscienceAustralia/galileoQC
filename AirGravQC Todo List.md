# AirGravQC Todo List

1. GRANTED - Get permissions from Xcalibur to use Canobie field data for tutorial.
2. Get permissions from SGL to use GA/GSV field data for tutorial.
3. Get permissions from Bell Geospace to use some field data for tutorial.
4. Get permissions from Andy Gabell to use some field data for tutorial.
5. Change all QC routines to report in the same format (Projectname, Blockname, summary [, detail text [, plot]]).
6. Standardise plots as much as possible (plot titles, axis labels, use of units).
7. Design geoWhizz v1.1, and then upgrade to it.
8. Test the various height checks: checkSafeClearance, checkClearance, checkVert; that they work as intended.
9. Is `checkDrape` useful or redundant?
10. `checkGNSS` only finds at most one error per line. It should find and report against the actual spec like `checkVertPlan`. Similarly it should have an option to plot lines with failures; and an option to test a specified array of lines.
11. All `check***` routines to have options `lines=[]`, `plot_flag=False`, `known=''`.
12. Standardise QC (check) functions to report: 1 - what it is doing; 2 - a one sentence result summary; 3 - line-by-line summary.
13. DONE. Stop using "erm" - use "grd" instead.


## geoWhizz v1.1

1. Get rid of separate CoordinateFrame and hold that info as project metadata.
2. Put in EPSG metadata.
3. Allow user to invent and store their own metadata at Project, Line, and Channel level. In other words, some metadata is there "just for the record" and is not used by the s/w and we can allow anything the user might want (and of course, `reportWhizz` will report it by design), and some metadata is actually used by the s/w and must follow the implicit conventions.
4. Stop storing channels that are expected, and are constant: `FLIGHT`, `LINE`, `PROJECTNUM`, `DOY`, etc. Channels constant across a line should be Line metadata and those constant across the project should be Project metadata. In every case, check that they are constant before reducing to metadata; report the error if they are not and do not convert to metadata in that case.
5. Channels that should be integers (e.g. observed number of GNSS satellites) must be checked that they are integers by `xyzToHDF`, and stored as integers.
6. Convert dates (e.g. DOY, decimal year, etc) to ISO date using Python's `DateTime` using function `doyToISO8601`.

## Other

1. Vertical axis of second plot needs to be in `xxx xxx m` format.
2. Dittio for checkDrape, checkClearance.
3. Ditto for checkSpeed.
4. Ditto for checkRepeats.
5. checkRepeats to be re-factored with the plot function (called four times) separate.
6. checkSpeed plots need horiz axis to match criteria (maxDistance -> easting or northing; maxDuration -> time, maxCount -> Fiducial).
7. DONE. All QC routines plotting a box-whisker plot to do it via a call to the same wp.plotBoxWhisker() function.
8. DONE. plotBoxWhisker to have a label for the x axis.
9. DONE. diffNoiseVturb to have a label_lines Bool parameter which labels each point on the plot so that one can identify outliers.
10. Ditto ilsNoiseVturb.
11. DONE. Have `allChanStats` capable of reporting statistics on channels after mean removal and/or first differencing.
12. checkSafeClearance needs a plot_flag parameter.
13. If there are lots of lines, and an error persists over all (or most) lines, then the report can be very long. Provide an `verbose=False` parameter which, if true, limits the report to the summary.
14. Calculating the statistics of angular variables has problems with 0, 2pi. Find a fix and implement.
15. `checkRawAGG` references the Lawin QC report. Replace this with an actual description in the Methods documentation.
