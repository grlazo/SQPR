#!/bin/perl

# -----------------------------------------------------------------------
#  File Name: runphred.cgi v.8                                 
#                                                       
#  Description: This Perl cgi-bin script reads from the 
#               Trace File Data Directory (TFDD) and 
#               executes the phred program. The output 
#               uses the directory name and appends
#		the output with *.fasta (FASTA sequence file), 
#		*.fasta.qual (FASTA quality file) and 
#               *.fasta.hist (phred histogram).
#                 
#  Created: 03/07/2000
#  Modified: 04/11/2000
# -----------------------------------------------------------------------

require "GlobalVariables.lib";
require "subparseform.lib";
&Parse_Form;

print "Content-type: text/html\n\n";

$dir = $formdata{'dir'};
$script = $ENV{'SCRIPT_NAME'};
$browser = $ENV{'HTTP_USER_AGENT'};

$fred = $PHRED_PARAMETER_FILE;
$fredenv = "$fred; export PHRED_PARAMETER_FILE";

$fredloc = "$PHRED";
$goseqs = "$TRACE_DIR";
$goseqd	= "$FASTA_DIR";

$fasta = "$goseqd/$dir.fasta";
$qual  = "$goseqd/$dir.fasta.qual";
$hist  = "$goseqd/$dir.fasta.hist";

$fastaext = "-sa $fasta";
$qualext  = "-qa $qual";
$histext  = "-qr $hist";

system "$fredenv; $fredloc -id $goseqs/$dir -trim_alt \"\"  $fastaext $qualext $histext&";
#system "sleep 90 |  chmod 666 $fasta $qual $hist&";

#----- HTML CODING--------

print <<"HEADER code";
<HTML>
<HEAD>
<TITLE>SQPR: Running phred analysis on $dir </TITLE>
</HEAD>

<BODY BGCOLOR="FFFFCC"> 
HEADER code

print <<"WAIT";
<CENTER>
<H2>SQPR: Running PHRED Analysis on <P><FONT COLOR="green"><I>$dir</I></FONT></H2><HR>
<H2>Please Wait! </H2>
<H3>This page will refresh shortly (depending on time require to run all trace samples).<BR>
    In the mean time, take a short break and relax. It should be less than a minute.<BR>
    Thank you.</H3>
</CENTER>
WAIT

if ($browser =~ /MSIE/) {
   print "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"40;URL=$CGI_BIN/fastadisplay.cgi?dir=$dir&script=$script\">";
 } else {
   print "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"5;URL=$CGI_BIN/fastadisplay.cgi?dir=$dir&script=$script\">";
 }

print <<"END";
<p><hr>

</BODY>
</HTML>
END

# -----------------------------------------------------------------------
# For more information see:
# Lazo, G.R., Tong, J., Miller, R., Hsia, C., Rausch, C., Kang, Y., 
# and Anderson, O.D. 2001. Software scripts for quality checking of 
# high-throughput nucleic acid sequencers. BioTechniques Vol 30, No. 6.
#
# This software is a "United States Government Work" under the terms 
# of the United States Copyright Act. It was written as part of official 
# duties performed by a United States Government employee(s) and thus 
# cannot be copyrighted. This software is freely available to the public 
# for use. The United States Department of Agriculture and the U.S. 
# Government have not placed any restriction on its use or reproduction. 
# -----------------------------------------------------------------------

