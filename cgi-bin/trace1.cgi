#!/bin/perl

# -----------------------------------------------------------------------
#  Filename: trace1.cgi v.8
#
#  Description: This Perl CGI script reads in the directory 
#		name and the trace file name.  With that 
#		trace file data, it will create a viewable 
#		trace file for the 384-format plate.
#
#  Created: 07/18/2000
#  Modified: 07/19/2000
# -----------------------------------------------------------------------

require "GlobalVariables.lib";
require "subparseform.lib"; 
&Parse_Form;

print "Content-type: text/html\n\n";

$directory = $formdata{'filename'};
$well = $formdata{'well'};
$start = $formdata{'start'};
#print "well = $well<BR>";
#print "start = $start<BR>";

$dir = "$TRACE_DIR/$directory";

$welloc = &deWELL($well,$start);			# Decode well location
#print "directory = $directory<BR>";
#print "welloc = $welloc<BR>";

$tracefile = &findtrace($directory,$welloc);		# find trace file with welloc
#print "trace = $tracefile<BR>";

$fname = "$WEB_TRACE_LINK_DIR/$tracefile";

system "rm -f $WEB_TRACE_DIR/*";		# Clear all the link in the www dir
system "ln -s $dir/$tracefile $WEB_TRACE_DIR/$tracefile";	#create a link to a www dir

print <<"HTML code";
<HTML>
<HEAD>
<TITLE> $tracefile </TITLE>
</HEAD>
<BODY BGCOLOR="FFFFCC"> 
<CENTER> $tracefile
<P>

<applet codebase="$JAVA_APPLET_CLASS_URL" code="ChromatogramApplet.class" width=400 height=200>
<param name="file" value="$fname">
</applet>

</CENTER>
</BODY>
</HTML>
HTML code

#---------Subroutines---------------------------
sub deWELL {
        $w = $_[0];		# well index
	$s = $_[1];		# starting well location for 384-format
	
	if (($w%16) == 1) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'A'.$num;
	} elsif (($w%16) == 2) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'B'.$num;
	} elsif (($w%16) == 3) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'C'.$num;
	} elsif (($w%16) == 4) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'D'.$num;
	} elsif (($w%16) == 5) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'E'.$num;
	} elsif (($w%16) == 6) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'F'.$num;
	} elsif (($w%16) == 7) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'G'.$num;
	} elsif (($w%16) == 8) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'H'.$num;
	} elsif (($w%16) == 9) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'I'.$num;
	} elsif (($w%16) == 10) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'J'.$num;
	} elsif (($w%16) == 11) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'K'.$num;
	} elsif (($w%16) == 12) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'L'.$num;
	} elsif (($w%16) == 13) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'M'.$num;
	} elsif (($w%16) == 14) {
		$num = sprintf("%02d", $w/16);
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'N'.$num;
	} elsif (($w%16) == 15) {
		$num = sprintf("%02d", $w/16);		
		$num = $num + $s;
		if ($num <= 9) {$num = '0'.$num;}
		return $welloc = 'O'.$num;
	} elsif (($w%16) == 0) {
		$num = sprintf("%02d", $w/16);
#		print "num = $num<BR>";
		$num = $num + $s - 1;
#		print "$num= $s<BR>";
		if ($num <= 9) {$num = '0'.$num;}
		$welloc = 'P'.$num;
#		print "welloc2 = $welloc<BR>";
		return $welloc;
	}
}

sub findtrace {
	$d = $_[0];		# directory name
	$w = $_[1];		# well location

        chdir ("$TRACE_DIR/$d");       
        
        opendir(TR, ".") || &ErrorMessage;
        @trfiles = readdir(TR);
        closedir(TR);
        shift(@trfiles); shift(@trfiles);
        if (@trfiles) {
           foreach $eatr (@trfiles) {
		$trace = $eatr;
                $eatr =~ tr/a-z/A-Z/;
                @nametr = split(/\_/, $eatr);
                $trwell = $nametr[1];
                if ($w eq $trwell) { return $tracefile = $trace;}        # get the trace file name
           }
        }   
}

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

