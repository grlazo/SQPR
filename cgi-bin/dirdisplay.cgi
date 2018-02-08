#!/bin/perl

# -----------------------------------------------------------------------
#  File Name: dirdisplay.cgi v.8                              
#                                                       
#  Description: This Perl cgi-bin script reads from the 
#               Trace File Data Directory (TFDD) and 
#               displays and links these directories. 
#               Selecting the directory initiates phred
#               processing with runphred.cgi.     
#                                                       
#  Created: 03/07/2000
#  Modified: 06/09/2000 modified style of the href
# -----------------------------------------------------------------------

require "GlobalVariables.lib";

$exist = 0;

chdir ("$TRACE_DIR");
opendir (DAT, ".") || &ErrorMessage;
@dirlist = readdir (DAT);
closedir (DAT);

chdir ("$FASTA_DIR");
opendir (DIR, ".") || &ErrorMessage;
@files = readdir (DIR);
closedir (DIR);

print "Content-type: text/html\n\n";

print <<"HEADER HTML";
<HTML>
<HEAD>
<TITLE>SQPR: Directories Containing Traces Ready for <b><i>phred</i></b> Processing</TITLE>
</HEAD>

<BODY BGCOLOR="FFFFCC">
<H1>SQPR:</H1>
<H2>Directories Containing Traces Ready for <b><i>phred</i></b> Processing</H2>
<HR>
HEADER HTML

shift(@files); shift(@files);
shift(@dirlist); shift(@dirlist);
@files = sort(@files);
@dirlist = sort(@dirlist);

$j = 0;

foreach (@files) {
 @name = split(/\./, $_);
 $filename[$j] = $name[0]; 
 $j++;
}
$numfile = $filename;

if (@dirlist) {
   print "<P><H3>Please Select A Directory:</H3>";
   foreach $dirt (@dirlist) {
      $exist = 0;
      @dname = split(/\./, $dirt);
      for ($i=0; $i <= $numfile; $i++) {
	 if ($dname[0] eq $filename[$i]) {
	     $exist = 1;
	 }
      }
      if ($exist == 0) {
	    print "<LI><A HREF=\"$CGI_BIN/runphred.cgi?dir=$dirt\" style=\"text-decoration:none\">$dirt</A><BR>\n" unless ($dirt =~ /\./);
      }
   }
}

print <<"FOOTER HTML";
</BODY>
</HTML>
FOOTER HTML

sub ErrorMessage {
	print "Content-type: text/html\n\n";
	print "The server can't access the file. Aborting script. \n";
	exit;
}

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
