#!/bin/perl

# -----------------------------------------------------------------------
#  File Name: phddisplay.cgi v.8                                  
#                                                       
#  Description: This Perl cgi-bin script reads from 
#               a selected *.qual in the Phred Output 
#		Data Directory (PODD) and displays 
#		analysis the file. The display is a 
#		histogram of Quality Read Length (QRL), 
#		a collection of sample measurements,  
#		and HTML color-coded tables for each 
#		sample for QRL, Phred Score, and more 
#		determined by data in the *.qual file
#                 
#  Created: 02/18/2000
#  Modified: 07/24/2000 v.6: added 384a & 384b plate display
#  Modified: 08/01/2000 v.6: modified histogram display
#  Modified: 09/25/2000 v.7: optimization
# -----------------------------------------------------------------------

require "GlobalVariables.lib";
require "subparseform.lib";
&Parse_Form;

print "Content-type: text/html\n\n";

$pass = $in = 0;
$file = $formdata{'file'};				# set the seq file to be read
$pass = $formdata{'pass'};				# set if seq file is pass or not

$qualfile = $file . '.qual';				# seq's qual file
$histfile = $file . '.hist';				# seq's hist file -- no longer used

$file1 = "$FASTA_DIR/$qualfile";			# path of the qual file
$file2 = "$FASTA_DIR/$histfile";			# path of the hist file

$passfile = "$FASTA_DIR/passfile.txt";
$failfile = "$FASTA_DIR/failfile.txt";

$avgstC = $avgqbC = $totalstC = $totalqbC = $qbC = $stC = 0;		# initialize variables
$avgstT = $avgqbT = $totalstT = $totalqbT = $qbT = $stT = 0;
$total20 = $totalavg = $avg20 = $avgAvg = $numbad = $numavg = $numzero = $numwell = 0;
$format96 = $format384a = $format384b = 0;
$i = $j = $ABI = $nodata = 0;
@wellC = ();

@set1 = ("A", "C", "E", "G", "I", "K", "M", "O");
@set2 = ("B", "D", "F", "H", "J", "L", "N", "P");
@set96 = ("A", "B", "C", "D", "E", "F", "G", "H");


system "chmod 666 $FASTA_DIR/$file $file1 $file2";	# change permission of the files
							# so they can be deleted if needed
@Sname = split(/\./, $file);
$Seqname = $Sname[0];

# Reading from Data File
open(INFO, $file1);
@lines = <INFO>;
close(INFO);

chomp @lines;						# remove all the carriage returns
$allSeq = "@lines";					# put all the element into 1 variable for parsing
@Seqbreak = split(/\>/, $allSeq);			# parse out each sequence
shift(@Seqbreak);					# remove the junk stuff in the beginning


foreach $eachSeq (@Seqbreak) {
   $totalbase = $badbase = 0;				# initialize the variables for each seq
   
   $numwell++;						# count number of sequences read
   @Seqparts = split(/ABI|SCF/, $eachSeq);		# split the seq header and the base value
   @Seqinfo = split(/\s+/, $Seqparts[0]);		# split the seq header into seq name, high, low, good
   @Seqname = split(/\_/, $Seqinfo[0]);			# split the seq name into well location, capillary value..etc
   
   $high = $Seqinfo[1];					# total of quality base for current seq
   $low = $Seqinfo[2];					# starting of good quality base for current seq
   $goood = $Qtotal = $Seqinfo[3];			# number of good quality bases
   	
   $well = $Seqname[1];					# storing well location
   $platename = $Seqname[0];				# storing the plate name
   
   $lastseqinfo = pop(@Seqname);			# getting the last set of info
   @lastinfo = split(/\./, $lastseqinfo);		# parsing out the capillary from cap.ab1
   $capname = $lastinfo[0];
   $ext = $lastinfo[1];			    		# getting the sequence extension

   if ($ext eq 'ab1') {$ABI = 1;}			# flag if ABI file
   elsif ($ext eq 'scf') {$ABI = 0;}			# flag if CEQ file
   
   $quall = $goood/($high-$low);			# calculate the quality index
   
   $qual{$well} = $quall;				# storing the sample quality index
   $fullname{$well} = $Seqinfo[0];			# storing the sequence name
   $cap{$well} = $capname;				# storing the capillary value
   $good{$well} = $goood;				# storing the good value
   
   if (($low == 0) && ($goood == 0)) {			# No-Good sample detected
      $bad = $badvalue{$well} = 100;			# seq is 100% bad
      $avg = $avgvalue{$well} = 0;			# none are good
      $numzero++;					# counting number of samples that are zero
   } else {						# Good sample detected
      @section = split(/\s+/, $Seqparts[1]);		# separate each base value
      for($i=$low+1; $i<=($low+$goood); ++$i) {
         $totalbase += $section[$i];			# adding all the good base together
         if ($section[$i] < 20) {$badbase++;}		# sum all the bad bases
      }
      $bad = (sprintf("%.2f", ($badbase/$Qtotal))*100);
      $avg = sprintf("%.0f", ($totalbase/$Qtotal));
      $badvalue{$well} = $bad;				# store bad value
      $avgvalue{$well} = $avg;				# store avg value
   }
   $totalavg += $totalbase;				# total value of all bases
   $total20 += $badbase;				# total value of bad bases
   $numavg += $Qtotal;					# total number of bases
   $numbad += $Qtotal;					# total number of bad bases

   if ($Seqinfo[0] =~ /\_...[zZ]/) {			# regular sample is detected
      $control = 0;
      if ($goood != 0) {
         $totalstT += $low;
         $totalqbT += $goood;
         $qbT++;
         $stT++;
      }
   } elsif ($Seqinfo[0] =~ /\_...[yY]/) {		# control sample is detected
      $control = 1;
      push(@wellC, $well);				# keep record of all the control well location
      if ($goood != 0) {
         $totalstC += $low;
         $totalqbC += $goood;
         $qbC++;
         $stC++;
      }
   }
   
   if (($well eq 'A01') || ($well eq 'a01')) {$A01 = 1;}	# flag if any target for plate format 
   elsif (($well eq 'A02') || ($well eq 'a02')) {$A02 = 1;}	# is detected
   elsif (($well eq 'A03') || ($well eq 'a03')) {$A03 = 1;}
   elsif (($well eq 'A04') || ($well eq 'a04')) {$A04 = 1;}
   elsif (($well eq 'A07') || ($well eq 'a07')) {$A07 = 1;}
   elsif (($well eq 'A08') || ($well eq 'a08')) {$A08 = 1;}
   elsif (($well eq 'A13') || ($well eq 'a13')) {$A13 = 1;}
   elsif (($well eq 'A14') || ($well eq 'a14')) {$A14 = 1;}
   elsif (($well eq 'A19') || ($well eq 'a19')) {$A19 = 1;}
   elsif (($well eq 'A20') || ($well eq 'a20')) {$A20 = 1;}
   elsif (($well eq 'A21') || ($well eq 'a21')) {$A21 = 1;}
   elsif (($well eq 'A22') || ($well eq 'a22')) {$A22 = 1;}
   elsif (($well eq 'A23') || ($well eq 'a23')) {$A23 = 1;}
   elsif (($well eq 'A24') || ($well eq 'a24')) {$A24 = 1;}
   elsif (($well eq 'B01') || ($well eq 'b01')) {$B01 = 1;}
   elsif (($well eq 'B02') || ($well eq 'b02')) {$B02 = 1;}
   elsif (($well eq 'B03') || ($well eq 'b03')) {$B03 = 1;}
   elsif (($well eq 'B04') || ($well eq 'b04')) {$B04 = 1;}
   elsif (($well eq 'B21') || ($well eq 'b21')) {$B21 = 1;}
   elsif (($well eq 'B22') || ($well eq 'b22')) {$B22 = 1;}
   elsif (($well eq 'B23') || ($well eq 'b23')) {$B23 = 1;}
   elsif (($well eq 'B24') || ($well eq 'b24')) {$B24 = 1;}
   elsif (($well eq 'I01') || ($well eq 'i01')) {$I01 = 1;}
   elsif (($well eq 'I02') || ($well eq 'i02')) {$I02 = 1;}
   elsif (($well eq 'I07') || ($well eq 'i07')) {$I07 = 1;}
   elsif (($well eq 'I08') || ($well eq 'i08')) {$I08 = 1;}
   elsif (($well eq 'I13') || ($well eq 'i13')) {$I13 = 1;}
   elsif (($well eq 'I14') || ($well eq 'i14')) {$I14 = 1;}
   elsif (($well eq 'I19') || ($well eq 'i19')) {$I19 = 1;}
   elsif (($well eq 'I20') || ($well eq 'i20')) {$I20 = 1;}
    	
}

if (($stC == 0) && ($qbC == 0)) {$avgstC = 'NA'; $avgqbC = 'NA';}  # no control sample
else {
   $avgstC = sprintf("%.0f", $totalstC/$stC);		# calculate the avg start position for control
   $avgqbC = sprintf("%.0f", $totalqbC/$qbC);		# calculate the avg read length for control
}

if (($stT == 0) && ($qbT == 0)) {$avgstT = 'NA'; $avgqbT = 'NA';}  # no test sample
else {
   $avgstT = sprintf("%.0f", $totalstT/$stT);		# calculate the avg start position for test sample
   $avgqbT = sprintf("%.0f", $totalqbT/$qbT);		# calculate the avg read length for test sample
}

if ($numavg == 0) {$avgAvg = 0; $avg20 = 0;}		# zero average
else {
   $avgAvg = (sprintf("%.0f", ($totalavg/$numavg)));	# calculate avg of all the avg value
   $avg20 = (sprintf("%.2f", ($total20/$numbad))*100);	# calculate avg of all the bad value
}

@allwell = keys(%fullname);				# getting a list of well name from the %fullname
@alwell = sort(@allwell);				# sort the well name
$first = shift(@alwell);				# getting the starting well location
$last = pop(@alwell);					# getting the ending well location
$startname = $fullname{$first};				# getting the name of the starting well location
$endname = $fullname{$last};				# getting the name of the ending well location

if (($startname) && ($endname)) {$nodata = 0;}		# there is data detected
elsif ((!$startname) && (!$endname)) {$nodata = 1;} 	# there is no data detected

if (($A01 == 1) && ($A02 == 1) && ($A07 == 1) && ($A08 == 1)) {$format96 = 1;} 		# 96 plate format detected
elsif (($A01 == 1) && ($A02 == 1) && ($I01 == 1) && ($I02 == 1)) {$format384a = 1;}	# 384-6x16 format detected
elsif (($A07 == 1) && ($A08 == 1) && ($I07 == 1) && ($I08 == 1)) {$format384a = 1;}
elsif (($A13 == 1) && ($A14 == 1) && ($I13 == 1) && ($I14 == 1)) {$format384a = 1;}
elsif (($A19 == 1) && ($A20 == 1) && ($I19 == 1) && ($I20 == 1)) {$format384a = 1;}
elsif (($A01 == 1) && ($A03 == 1) && ($A21 == 1) && ($A23 == 1)) {$format384b = 1;}	# 384-8x12 format detected
elsif (($A02 == 1) && ($A04 == 1) && ($A22 == 1) && ($A24 == 1)) {$format384b = 1;}
elsif (($B01 == 1) && ($B03 == 1) && ($B21 == 1) && ($B23 == 1)) {$format384b = 1;}
elsif (($B02 == 1) && ($B04 == 1) && ($B22 == 1) && ($B24 == 1)) {$format384b = 1;}

foreach (@allwell) {
   $flag = 0;
   
   @nameL = split(/\d/, $_);
   @nameD = split(/\D/, $_);
   $wL = $nameL[0];				# well letter
   $wD = $nameD[1];				# well number

   if ($format96 == 1) {			# assign index num in 96-format
      $wD = &defrow($wL,$wD);
   } elsif ($format384b == 1) {			# assign index in 384b-format
      $wD = &def384brow($wL,$wD);
   } else {					# assign index num in 384-format
      $wD = &def384row($wL, $wD);
   }
#   $ind = $db{$_};				# getting percentage
   $ind = $qual{$_};				# getting percentage

   chomp($_);
   foreach $wellCtrl (@wellC) {			# check if well is control or not
      if ($wellCtrl eq $_) {			# flag if it is a control
	 $flag = 1;
      }
   }
   if ($flag == 1) {
       $wC[$wD] = 1;				# set to 1 if control
   } else {
       $wC[$wD] = 0;				# set to 0 if not control 
   }
   $val = &defcolor($ind);

   $number[$wD] = (sprintf("%.2f", $ind) * 100);
   $color[$wD] = $val;				# set color
}

#----Capillary Value Array-----
if ($ABI == 1) {
 while (($wel, $capvalue) = each (%cap)) {
   @capL = split(/\d/, $wel);
   @capD = split(/\D/, $wel);
   $cL = $capL[0];
   $cD = $capD[1];

   if ($format96 == 1) {			# assign well index in 96-format
      $n = &defrow($cL,$cD);
   } elsif ($format384b == 1) {			# assign well index in 384b-format
      $n = &def384brow($cL,$cD);
   } else {
      $n = &def384row($cL,$cD);
   }
   $Capillary[$n] = $capvalue;			# store capillary value
 }
}

#-----Good Base Value Array-----
while (($wel, $govalue) = each (%good)) {
   @goL = split(/\d/, $wel);
   @goD = split(/\D/, $wel);
   $gL = $goL[0];
   $gD = $goD[1];

   if ($format96 == 1) {			# assign well index in 96-format
      $z = &defrow($gL,$gD);
   } elsif ($format384b == 1) {			# assign well index in 384b-format
      $z = &def384brow($gL,$gD);
   } else {
      $z = &def384row($gL,$gD);
   }	
#   print "z = $z ; goodvalue = $govalue<BR>";		
   $Goodvalue[$z] = $govalue;			# store good value
   $Goodcolor[$z] = &deGcolour($govalue);	# storing the bg color of cell
}						

#-----Bad Base Value Array-----
while (($wel, $bdvalue) = each (%badvalue)) {
   @goL = split(/\d/, $wel);
   @goD = split(/\D/, $wel);
   $gL = $goL[0];
   $gD = $goD[1];

   if ($format96 == 1) {			# assign well index in 96-format
      $z = &defrow($gL,$gD);
   } elsif ($format384b == 1) {			# assign well index in 384b-format
      $z = &def384brow($gL,$gD);
   } else {
      $z = &def384row($gL,$gD);
   }	
#   print "z = $z ; badvalue = $bdvalue<BR>";		
   $badvalue[$z] = $bdvalue;			# store good value
   $badcolor[$z] = &defbadcolor($bdvalue);	# storing the bg color of cell
}					

#-----Avg Base Value Array-----
while (($wel, $avvalue) = each (%avgvalue)) {
   @goL = split(/\d/, $wel);
   @goD = split(/\D/, $wel);
   $gL = $goL[0];
   $gD = $goD[1];

   if ($format96 == 1) {			# assign well index in 96-format
      $z = &defrow($gL,$gD);
   } elsif ($format384b == 1) {			# assign well index in 384b-format
      $z = &def384brow($gL,$gD);
   } else {
      $z = &def384row($gL,$gD);
   }	
#   print "z = $z ; avgvalue = $avvalue<BR>";		
   $avgvalue[$z] = $avvalue;			# store good value
   $avgcolor[$z] = &defavgcolor($avvalue);	# storing the bg color of cell
}					

#----- Hist.cgi-------------

$highest = $total300 = $total350 = 0;
$xhigh = 15;
@GV = @Goodvalue;
shift(@GV);

foreach (@GV) {
   if (($_ >= 0) && ($_ <= 49))		{$hist[0]++;}			# count values between 0-49
   if (($_ >= 50) && ($_ <= 99))	{$hist[1]++;}					# between 50-99
   if (($_ >= 100) && ($_ <= 149))	{$hist[2]++;}					# between 100-149
   if (($_ >= 150) && ($_ <= 199))	{$hist[3]++;}					# between 150-199
   if (($_ >= 200) && ($_ <= 249))	{$hist[4]++;}					# between 200-249
   if (($_ >= 250) && ($_ <= 299))	{$hist[5]++;}					# between 250-299
   if (($_ >= 300) && ($_ <= 349))	{$hist[6]++; $total300 ++;}			# between 300-349
   if (($_ >= 350) && ($_ <= 399))	{$hist[7]++; $total300 ++; $total350 ++;}	# between 350-399
   if (($_ >= 400) && ($_ <= 449))	{$hist[8]++; $total300 ++; $total350 ++;}	# between 400-449
   if (($_ >= 450) && ($_ <= 499))	{$hist[9]++; $total300 ++; $total350 ++;}	# between 450-499
   if (($_ >= 500) && ($_ <= 549))	{$hist[10]++; $total300 ++; $total350 ++;}	# between 500-549
   if (($_ >= 550) && ($_ <= 599))	{$hist[11]++; $total300 ++; $total350 ++;}	# between 550-599
   if (($_ >= 600) && ($_ <= 649))	{$hist[12]++; $total300 ++; $total350 ++;}	# between 600-649
   if (($_ >= 650) && ($_ <= 699))	{$hist[13]++; $total300 ++; $total350 ++;}	# between 650-699
   if (($_ >= 700) && ($_ <= 749))	{$hist[14]++; $total300 ++; $total350 ++;}	# between 700-749
   if ($_ >= 750)			{$hist[15]++; $total300 ++; $total350 ++;}	# between 750+
}

foreach (@hist) {
   if ($highest < $_) {$highest = $_;} 			# find the highest count among the 50 groups
}


#-------Passing On?---------------------------
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);	# using the date function
$year += 1900;
$mon += 1;
$mon = sprintf("%02d", $mon);
$day = sprintf("%02d", $mday);						# to get the current date


$date = $year . '-' . $mon . '-' . $day;

if($pass eq 'Yes') {
   open(PS, $passfile);							# getting pass content
   @passcontent = <PS>;
   close(PS);

   foreach (@passcontent) {
	chomp($_);
	@filecontent = split(/\s/, $_);
	if ($Seqname eq $filecontent[1]) {
	   $in = 1;
 	}
   }
   if ($in == 0) {
  	open(PF, ">>$passfile");
   	print PF "$date $Seqname $platename $total300 $total350 $avgqbT $avgAvg $avg20 $numzero $numavg\n";
   	close(PF);
   }
} elsif($pass eq 'No') {	 
   open(FA, $failfile);
   @failcontent = <FA>;
   close(FA);

   foreach (@failcontent) {
   	chomp($_);
	@filecontent = split(/\s/, $_);
	if ($Seqname eq $filecontent[1]) {
	   $in = 1;
 	}
   }
   if ($in == 0) {
  	open(FF, ">>$failfile");
   	print FF "$date $Seqname $platename $total300 $total350 $avgqbT $avgAvg $avg20 $numzero $numavg\n";
   	close(FF);
   }
}

#------End of Passing On----------------------


#----- Header for the HTML Page---------------

print <<"HEADER code";
<HTML>
<HEAD>
<TITLE> $Seqname </TITLE>

<script LANGUAGE="Javascript1.2">
   function MessageScreen(file,tar,w,h) 
	{window.open(file,tar,"width="+w+",height="+h+"");}
</script>

</HEAD>
<BODY BGCOLOR="FFFFCC"> 
HEADER code
#-----End of Header for the HTML Page---------


#-------- Output Histigram to the screen-------
print <<"HIST TITLE";
<H4><U>Base Quality Histogram for:  <I>$Seqname</I> (Plate: $platename)</U></H4>

<TABLE BORDER=1 COLS=2>
<TR>
   <TD ROWSPAN=9>
      <FONT SIZE=-6>Y-axis: # of Samples</FONT><BR>
      <TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>
HIST TITLE

for ($i=$highest; $i>0; $i--) {
   print "<TR><TH WIDTH=50 HEIGHT=15>$i</TH>";
   for ($j=0; $j<=$xhigh; $j++) {
	if($hist[$j] >= $i) {
		print "<TD BGCOLOR=\"#000000\" HEIGHT=15 WIDTH=15 ALIGN=\"center\">*</TD>";
	}
	else {
		print "<TD HEIGHT=15 WIDTH=15>&nbsp;</TD>";
	}
	print "<TD HEIGHT=15 WIDTH=5>&nbsp;</TD>";
   }
   print "</TR>\n";
}

print "<TR><TD WIDTH=50 HEIGHT=15>&nbsp;</TD>";

$xvalue = 0;
for ($i=0; $i<$xhigh+1; $i++) {
   @len = split(//, $xvalue);
   print "<TD WIDTH=15 HEIGHT=15 VALIGN=\"top\"><B>";
   foreach (@len) {
	print "$_<BR>";
   }
   print "</B></TD><TD HEIGHT=15 WIDTH=5>&nbsp;</TD>";
   $xvalue += 50;
}

print "</TR>\n";
print "</TABLE>\n";
#------End of Output of Histogram--------------


#------OUTPUT OF HIST&OTHER INFO--------------
print <<"HIST&OTHER INFO";
      <FONT SIZE=-6>X-axis: # of Quality Bases</FONT>
   </TD>

   <TD WIDTH="100" COLSPAN=3 VALIGN="bottom">
   <FONT SIZE=-5>Note: &nbsp; Each block (shown left) represents a </FONT><BR>
   <FONT SIZE=-5> &nbsp;&nbsp; sample in a range of quality bases </FONT><BR>
   <FONT SIZE=-5> &nbsp;&nbsp; using the phred histogram (-qr) option.</FONT><P>
   <FONT SIZE=-5> The command line for generating files used in the following output was: </FONT><P>
   <FONT SIZE=-5> phred -id &nbsp;&nbsp; <b>trace file data directory</b> </FONT><BR>
   <FONT SIZE=-5> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -trim_alt &nbsp; ""</FONT><BR>
   <FONT SIZE=-5> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -sa &nbsp; <b>fasta sequence output</b></FONT><BR>
   <FONT SIZE=-5> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -qa &nbsp; <b>fasta quality output</b></FONT><BR>
   <FONT SIZE=-5> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -qr &nbsp; <b>read lengths histogram</b></FONT><P>
   <FONT SIZE=-5> See: Ewing, Hillier, Wendl, and Green.</FONT><BR>
   <FONT SIZE=-5> &nbsp;&nbsp; 1998. Genome Research 8:175-185.</FONT><P>

   <FONT SIZE=-5>Number of Samples <u>&gt;</u> 300 = </FONT>&nbsp;<FONT SIZE=3><B> $total300 </B></FONT><BR>
   <FONT SIZE=-5>Number of Samples <u>&gt;</u> 350 = </FONT>&nbsp;<FONT SIZE=3><B> $total350 </B></FONT><BR><BR>
   </TD>
</TR>

<TR>
   <TD><FONT SIZE=2><B>Quality Base<BR>Averages for:</FONT></B></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>Controls</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>Test Samples</FONT></TD>
</TR>
<TR> 
   <TD><FONT SIZE=2>Start Position</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>$avgstC</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>$avgstT</FONT></TD>
</TR>
<TR>
   <TD><FONT SIZE=2>Read Length</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>$avgqbC</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>$avgqbT</FONT></TD>
</TR>
<TR>
   <TD colspan=3 align=center><FONT SIZE=2><b>Other Statistics</b></FONT></TD>
</TR>
<TR>
   <TD colspan=2><FONT SIZE=2>Number of Failed Samples</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2><b>$numzero</b> / $numwell</FONT></TD>
</TR>
<TR>
   <TD colspan=2><FONT SIZE=2>Percent Calls &lt; <i>phred20</i></FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>$avg20</FONT></TD>
</TR>
<TR>
   <TD colspan=2><FONT SIZE=2>Average <i>phred</i> Score</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>$avgAvg</FONT></TD>
</TR>
<TR>
   <TD colspan=2><FONT SIZE=2>Total Bases Counted</FONT></TD>
   <TD ALIGN=CENTER><FONT SIZE=2>$numavg</FONT></TD>
</TR>
</TABLE>
<P>
HIST&OTHER INFO

# Avg Bad ($total20/$numbad) = $avg20<BR>
# Avg Avg ($totalavg/$numavg) = $avgAvg<BR>
# Zero Well = $numzero/$numwell<BR>

if ($in == 0) {
print <<"HIST&OTHER INFO";
<FORM ACTION="$CGI_BIN/phddisplay.cgi" METHOD=POST>
<B> Should these sequences be passed on for further processing?</B>
	<INPUT TYPE="hidden" name="filename" value="$Seqname">
	<INPUT TYPE="hidden" name="file" value="$file">
	<input type="submit" name="pass" value="Yes">
	<input type="submit" name="pass" value="No">
	<input type="submit" name="pass" value="Don't Know">
</FORM>
HIST&OTHER INFO
}


print <<"HIST&OTHER INFO";
<P>
<TABLE BORDER=0>
   <TR>	<TD><B> Sample Range </B></TD>
	<TD> is $first : </TD>
	<TD><I> $startname </I></TD>
   </TR>
   <TR>	<TD>&nbsp; </TD>
	<TD> to $last : </TD>
	<TD><I> $endname </I></TD>
   </TR>
</TABLE>

<FORM>
<B>Comments:</B><BR>
<TEXTAREA NAME="comment" COLS=60 ROWS=5 WRAP></TEXTAREA>
</FORM>
HIST&OTHER INFO
#------End of HIST&OTHER INFO--------------

#--------INCOMPLETE Data? DELETE it------------
if ($numwell < 96) {
print <<"REDO";
<FORM ACTION="redo.cgi?file=$file" METHOD=POST>
<P>The following plate appears incomplete; do you wish to <FONT COLOR="blue">RE-RUN</FONT> this later:
<INPUT TYPE="submit" VALUE="Proceed?" ALIGN="middle">
</FORM>
REDO
}
#--------End of INCOMPLETE -> DELETE FORM---------------------

if ($nodata == 1) {			# cannot read the qual data, different format

print <<"NODATA";
<BR>
<CENTER>
<HR>
<H4> Cannot Display the Well-tables. <BR> No Data has been read. </H4>
<FONT SIZE=3>
Your Sequence name may not be in the correct format. <BR>
Please check the <A HREF=\"/homepage/wheatbiotech/abiname.html\">ABI3700 Naming System</A> or the <A HREF=\"/homepage/wheatbiotech/ceqname.html\">CEQ2000 Naming System</A><BR>	
For further assistance, please contact Gerry or Jenny. <BR>
 Thank you. 
</FONT>
</CENTER>
NODATA

} elsif ($nodata == 0) {		# able to read the qual data

   if ($format96 == 1) {		# display for 96-wells plate

#------Print QUality Read Lengths Table--------
print <<"QROW1";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<H4><U>Quality Read Lengths (<b>QRL</b>s) for: <I>$Seqname</I></U></H4>
QROW1

print <<"QTABLE1";
<CENTER>
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
QTABLE1

   for($i=1; $i <= 12; $i++) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

$k = 0;
foreach (@set96) {
print <<"QTABLE";
</TR>
<TR>
   <TD><CENTER> $_ </CENTER></TD>
QTABLE

   for($i=1+$k; $i <= 12+$k; $i++) {
     $bcolor = $Goodcolor[$i];
     $num = $Goodvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   $k += 12;
}

print <<"QTABLE10";
</TR>
</TABLE>
</CENTER>
QTABLE10

print <<"QBEG ROW2";
</TD>
<TD ALIGN=CENTER WIDTH=100>
QBEG ROW2

print <<"QKEY TABLE";
<TABLE BORDER=1 COLOR="000000" COLS=1 WIDTH=90>
<TR><TD ALIGN=CENTER> Color Key </TD></TR>
<TR><TD BGCOLOR = #FF66FF ALIGN=CENTER> 0 - 150 </TD></TR> 
<TR><TD BGCOLOR = #FF99FF ALIGN=CENTER> 151 - 300 </TD></TR> 
<TR><TD BGCOLOR = #33CCFF ALIGN=CENTER> 301 - 450 </TD></TR>
<TR><TD BGCOLOR = #3333FF ALIGN=CENTER><FONT COLOR=WHITE> 451 - 600 </FONT></TD></TR> 
<TR><TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE> 601+ </FONT></TD></TR>
</FONT></TD></TR>
</TABLE>
QKEY TABLE

print <<"QEND ROW2";
</TD>
</TR>
</TABLE>
QEND ROW2
#------End of Print Quality Base Table----

#------Print Percentage Bases Below 20 Table----------
print <<"ROW";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<H4><U>Percent Base Calls Below <I>Phred20</I> for: <I>$Seqname</I></U></H4>
ROW

print <<"BAD TABLE";
<CENTER>
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
BAD TABLE

   for($i=1; $i <= 12; $i++) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>
   <TD><CENTER> A </CENTER></TD>
BAD TABLE

   for($i=1; $i <= 12; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>
   <TD><CENTER> B </CENTER></TD>
BAD TABLE

   for($i=13; $i <= 24; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>   
   <TD><CENTER> C </CENTER></TD>
BAD TABLE

   for($i=25; $i <= 36; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>   
   <TD><CENTER> D </CENTER></TD>
BAD TABLE

  for($i=37; $i <= 48; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>   
   <TD><CENTER> E </CENTER></TD>
BAD TABLE

  for($i=49; $i <= 60; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>   
   <TD><CENTER> F </CENTER></TD>
BAD TABLE

  for($i=61; $i <= 72; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>   
   <TD><CENTER> G </CENTER></TD>
BAD TABLE

  for($i=73; $i <= 84; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
<TR>   
   <TD><CENTER> H </CENTER></TD>
BAD TABLE

  for($i=85; $i <= 96; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }

print <<"BAD TABLE";
</TR>
</TABLE>
</CENTER>
BAD TABLE

print <<"KEY TABLE";
</TD>
<TD ALIGN=CENTER WIDTH=100>
<TABLE BORDER=1 COLOR="000000" COLS=1 WIDTH=90>
<TR><TD ALIGN=CENTER> Color Key </TD></TR>
<TR><TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE> 0% - 20%</FONT></TD></TR> 
<TR><TD BGCOLOR = #33CCFF ALIGN=CENTER> 21% - 30% </TD></TR>
<TR><TD BGCOLOR = #FF99FF ALIGN=CENTER> 31% - 40% </TD></TR>
<TR><TD BGCOLOR = #FF66FF ALIGN=CENTER> 41%+ </TD></TR>
</TABLE>
</TD>
</TR>
</TABLE>
KEY TABLE

#---------End of Percentage Bases below 20 Table---------------

#------Print Average Phred score in Qread-length Table----------
print <<"ROW";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<H4><U>Average Phred Scores in <b>QRL</b>s for: <I>$Seqname</I></U></H4>
ROW

print <<"AVG TABLE";
<CENTER>
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
AVG TABLE

   for($i=1; $i <= 12; $i++) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

$k = 0;
foreach (@set96) {
print <<"AVG TABLE";
</TR>
<TR>
   <TD><CENTER> $_ </CENTER></TD>
AVG TABLE

   for($i=1+$k; $i <= 12+$k; $i++) {
     $bcolor = $avgcolor[$i];
     $num = $avgvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   $k += 12;
}

print <<"AVG TABLE";
</TR>
</TABLE>
</CENTER>
AVG TABLE

print <<"KEY TABLE";
</TD>
<TD ALIGN=CENTER WIDTH=100>
<TABLE BORDER=1 COLOR="000000" COLS=1 WIDTH=90>
<TR><TD ALIGN=CENTER> Color Key </TD></TR>
<TR><TD BGCOLOR = #FF66FF ALIGN=CENTER> 0 - 20 </TD></TR> 
<TR><TD BGCOLOR = #33CCFF ALIGN=CENTER> 21 - 30 </TD></TR>
<TR><TD BGCOLOR = #3333FF ALIGN=CENTER><FONT COLOR=WHITE> 31 - 40 </FONT></TD></TR>
<TR><TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE> 41+ </FONT></TD></TR>
</TABLE>
</TD>
</TR>
</TABLE>
KEY TABLE

#---------End of Average Phred score in Qread-length Table---------------

#---------Print Capillary Value Table-----------
if ($ABI == 1) {
print <<"CAP TABLE";
<BR>
<H4><U>Capillaries Used for: <I>$Seqname</I></U></H4>
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
CAP TABLE

   for($i=1; $i <= 12; $i++) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

$k = 0;
foreach (@set96) {
print <<"CAP TABLE";
</TR>
<TR>
   <TD><CENTER> $_ </CENTER></TD>
CAP TABLE

   for($i=1+$k; $i <= 12+$k; $i++) {
     $bcolor = $avgcolor[$i];
     $num = $Capillary[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- underline & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
        
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><FONT COLOR=$fcolor>$format</FONT></TD>\n";
   }
   $k += 12;
}

print <<"CAP TABLE";
</TR>
</TABLE>
CAP TABLE
}		# end of ABI-Capillary value
#-----------End of Capillary Value Table--------

  } elsif ($format384b == 1) {			# display for 384b-format (8x12)
  
  @str = split(/\_/, $startname);		# extract the well pos
  @sD = split(/\D/, $str[1]);			# extrat starting number pos
  @sL = split(/\d/, $str[1]);			# extrat starting letter pos
  $startn = $sD[1];				# store the starting well number
  $startn = $startn+1-1;			#    change character to integer
  $startl = $sL[0];				# store the starting well letter
  
  if (($startl eq 'A') || ($startl eq 'a')) {@setletter = @set1;}
  if (($startl eq 'B') || ($startl eq 'b')) {@setletter = @set2;}
     
  #------Print QUality Read Lengths Table in 8x12 384 format--------
print <<"Q384B";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<H4><U>Quality Read Lengths (<b>QRL</b>s) for: <I>$Seqname</I></U></H4>
Q384B

print <<"Q384B";
<CENTER>
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
Q384B

   for($i=$startn; $i <= 24; $i=$i+2) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

print <<"Q384B";
</TR>
Q384B

$j = 1;
foreach (@setletter) {
   print "<TR>\n";
   print "   <TD><CENTER> $_ </CENTER></TD>\n";

   for($i=$j; $i <= 12+($j-1); $i++) {
     $bcolor = $Goodcolor[$i];
     $num = $Goodvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace2.cgi?filename=$Seqname&well=$i&startn=$startn&startl=$startl";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   print"</TR>\n";
   $j = $i; 
}

print <<"Q384B";
</TABLE>
</CENTER>
Q384B

print <<"QBEG ROW2";
</TD>
<TD ALIGN=CENTER WIDTH=100>
QBEG ROW2

print <<"QKEY TABLE";
<TABLE BORDER=1 COLOR="000000" COLS=1 WIDTH=90>
<TR><TD ALIGN=CENTER> Color Key </TD></TR>
<TR><TD BGCOLOR = #FF66FF ALIGN=CENTER> 0 - 150 </TD></TR> 
<TR><TD BGCOLOR = #FF99FF ALIGN=CENTER> 151 - 300 </TD></TR> 
<TR><TD BGCOLOR = #33CCFF ALIGN=CENTER> 301 - 450 </TD></TR>
<TR><TD BGCOLOR = #3333FF ALIGN=CENTER><FONT COLOR=WHITE> 451 - 600 </FONT></TD></TR> 
<TR><TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE> 601+ </FONT></TD></TR>
</FONT></TD></TR>
</TABLE>
QKEY TABLE

print <<"QEND ROW2";
</TD>
</TR>
</TABLE>
QEND ROW2
#------End of Print Quality Base Table in 384b-format (8x12)----

#------Print Percentage Bases Below 20 Table in 384b format (8x12)----------
print <<"ROW";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<H4><U>Percent Base Calls Below <I>Phred20</I> for: <I>$Seqname</I></U></H4>
<P>
ROW

print <<"B384B";
<CENTER>
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
B384B

   for($i=$startn; $i <= 24; $i=$i+2) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

print <<"B384B";
</TR>
B384B

$j = 1;
foreach (@setletter) {
   print "<TR>\n";
   print "   <TD><CENTER> $_ </CENTER></TD>\n";

   for($i=$j; $i <= 12+($j-1); $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   print"</TR>\n";
   $j = $i; 
}

print <<"B384B";
</TABLE>
</CENTER>
B384B

print <<"KEY TABLE";
</TD>
<TD ALIGN=LEFT WIDTH=400>
<TABLE BORDER=1 COLOR="000000" COLS=1 WIDTH=90>
<TR><TD ALIGN=CENTER><FONT SIZE=2> Color Key </FONT></TD></TR>
<TR><TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE> 0% - 20% </FONT></TD></TR>
<TR><TD BGCOLOR = #33CCFF ALIGN=CENTER> 21% - 30% </TD></TR>
<TR><TD BGCOLOR = #FF99FF ALIGN=CENTER> 31% - 40%</TD></TR>
<TR><TD BGCOLOR = #FF66FF ALIGN=CENTER> 41%+ </TD></TR>
</TABLE>
</TD>
</TR>
</TABLE>
KEY TABLE

#---------End of Percentage Bases below 20 Table in 384b format (8x12)---------------

#------Print Average Phred score in Qread-length Table in 384b format (8x12)----------
print <<"A384B";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<H4><U>Average Phred Scores in <b>QRL</b>s for: <I>$Seqname</I></U></H4>
A384B

print <<"A384B";
<CENTER>
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
A384B

   for($i=$startn; $i <= 24; $i=$i+2) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

print <<"A384B";
</TR>
A384B

$j = 1;
foreach (@setletter) {
   print "<TR>\n";
   print "   <TD><CENTER> $_ </CENTER></TD>\n";

   for($i=$j; $i <= 12+($j-1); $i++) {
     $bcolor = $avgcolor[$i];
     $num = $avgvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace.cgi?filename=$Seqname&well=$i";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   print"</TR>\n";
   $j = $i; 
}

print <<"A384B";
</TABLE>
</CENTER>
A384B

print <<"KEY TABLE";
</TD>
<TD ALIGN=CENTER WIDTH=100>
<TABLE BORDER=1 COLOR="000000" COLS=1 WIDTH=90>
<TR><TD ALIGN=CENTER> Color Key </TD></TR>
<TR><TD BGCOLOR = #FF66FF ALIGN=CENTER> 0 - 20 </TD></TR> 
<TR><TD BGCOLOR = #33CCFF ALIGN=CENTER> 21 - 30 </TD></TR>
<TR><TD BGCOLOR = #3333FF ALIGN=CENTER><FONT COLOR=WHITE> 31% - 40% </FONT></TD></TR>
<TR><TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE> 41+ </FONT></TD></TR>
</TABLE>
</TD>
</TR>
</TABLE>
KEY TABLE

#---------End of Average Phred score in Qread-length Table in 384b format (8x12)---------------

#------Print Capillary score in Qread-length Table in 384b format (8x12)----------
print <<"C384B";
<BR>
<H4><U>Capillaries Used for: <I>$Seqname</I></U></H4>
C384B

print <<"C384B";
<TABLE BORDER=1 COLOR="000000" COLS=13 WIDTH=500>
<TR>   
   <TD>&nbsp;</TD>
C384B

   for($i=$startn; $i <= 24; $i=$i+2) {
	print "<TD><CENTER><FONT COLOR=BLACK> $i </FONT></CENTER></TD>\n";
   }

print <<"C384B";
</TR>
C384B

$j = 1;
foreach (@setletter) {
   print "<TR>\n";
   print "   <TD><CENTER> $_ </CENTER></TD>\n";

   for($i=$j; $i <= 12+($j-1); $i++) {
     $bcolor = $avgcolor[$i];
     $num = $Capillary[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- underline & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
        
     print "<TD ALIGN=CENTER BGCOLOR = $bcolor><FONT COLOR=$fcolor>$format</FONT></TD>\n";
   }
   print"</TR>\n";
   $j = $i; 
}

print <<"C384B";
</TABLE>
C384B
#---------End of Capillary score in Qread-length Table in 384b format (8x12)---------------

  } elsif ($format384a == 1) {			# display for 384-plate (6x16)
  
  @str = split(/\_/, $startname);		# extract the well pos
  @sD = split(/\D/, $str[1]);			# extrat starting number pos
  $start = $sD[1];				# store the starting well position
  $start = $start+1-1;				# change character to integer
  
#------Print QUality Read Lengths Table in 384 format--------
print <<"QROW1";
<BR>
<TABLE BORDER=0 COLS=1 WIDTH=600>
<TR>
<TD WIDTH=600>
<FONT SIZE=3><B><U>Quality Read Lengths (<b>QRL</b>s) for: <I>$Seqname</I></U></B></FONT><FONT SIZE=1> (Plate: $platename)</FONT>
<P>
QROW1

print <<"QTABLE";
<TABLE BORDER=1 COLOR="000000" COLS=17 WIDTH=600>
QTABLE

$k = 5; $j = 81;
while ($k >= 0) {
$Q = $start+$k; 
print <<"QTABLE";
<TR>
   <TD WIDTH=20><CENTER> $Q </CENTER></TD>
QTABLE

   for($i=$j; $i <= $j+15; $i++) {
     $bcolor = $Goodcolor[$i];
     $num = $Goodvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace1.cgi?filename=$Seqname&well=$i&start=$start";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR=$bcolor WIDTH=20><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   $k--;
   $j -= 16;
   print "</TR>\n";
}

print <<"QTABLE";
<TR>   
   <TD WIDTH=20>&nbsp;</TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> A </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> B </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> C </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> D </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> E </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> F </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> G </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> H </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> I </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> J </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> K </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> L </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> M </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> N </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> O </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> P </FONT></CENTER></TD>
</TR>
</TABLE>
QTABLE

print <<"QBEG ROW2";
</TD></TR>
<TR>
<TD ALIGN=LEFT WIDTH=600>
QBEG ROW2

print <<"QKEY TABLE";
<TABLE BORDER=1 COLOR="000000" COLS=6 WIDTH=400>
<TR><TD ALIGN=LEFT><FONT SIZE=2> Color Key </FONT></TD>
    <TD BGCOLOR = #FF66FF ALIGN=CENTER><FONT SIZE=2> 0 - 150 </FONT></TD>
    <TD BGCOLOR = #FF99FF ALIGN=CENTER><FONT SIZE=2> 151 - 300 </FONT></TD>
    <TD BGCOLOR = #33CCFF ALIGN=CENTER><FONT SIZE=2> 301 - 450 </FONT></TD>
    <TD BGCOLOR = #3333FF ALIGN=CENTER><FONT COLOR=WHITE SIZE=2> 451 - 600 </FONT></TD>
    <TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE SIZE=2> 601+ </FONT></TD></TR>
</TD></TR>
</TABLE>
QKEY TABLE

print <<"QEND ROW2";

</TD>
</TR>
</TABLE>
QEND ROW2
#------End of Print Quality Base Table in 384 format----

#------Print Percentage Bases Below 20 Table in 384 format----------
print <<"ROW";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<FONT SIZE=3><B><U>Percent Base Calls Below <I>Phred20</I> for: <I>$Seqname</I></U></B></FONT><FONT SIZE=1> (Plate: $platename)</FONT>
<P>
ROW

print <<"BAD TABLE";
<TABLE BORDER=1 COLOR="000000" COLS=17 WIDTH=600>
BAD TABLE

$k = 5; $j = 81;
while ($k >= 0) {
$B = $start+$k; 
print <<"BAD TABLE";
<TR>
   <TD WIDTH=20><CENTER> $B </CENTER></TD>
BAD TABLE

   for($i=$j; $i <= $j+15; $i++) {
     $bcolor = $badcolor[$i];
     $num = $badvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace1.cgi?filename=$Seqname&well=$i&start=$start";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR=$bcolor WIDTH=20><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   $k--;
   $j -= 16;
   print "</TR>\n";
}

print <<"BAD TABLE";
<TR>   
   <TD WIDTH=20>&nbsp;</TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> A </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> B </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> C </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> D </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> E </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> F </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> G </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> H </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> I </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> J </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> K </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> L </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> M </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> N </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> O </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> P </FONT></CENTER></TD>
</TR>
</TABLE>
BAD TABLE

print <<"KEY TABLE";
</TD></TR>
<TR>
<TD ALIGN=LEFT WIDTH=400>
<TABLE BORDER=1 COLOR="000000" COLS=5 WIDTH=400>
<TR><TD ALIGN=CENTER><FONT SIZE=2> Color Key </FONT></TD>
    <TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE SIZE=2> 0% - 20%</FONT></TD>
    <TD BGCOLOR = #33CCFF ALIGN=CENTER><FONT SIZE=2> 21% - 30% </FONT></TD>
    <TD BGCOLOR = #FF99FF ALIGN=CENTER><FONT SIZE=2> 31% - 40% </FONT></TD>
    <TD BGCOLOR = #FF66FF ALIGN=CENTER><FONT SIZE=2> 41%+ </FONT></TD></TR>
</TABLE>
</TD>
</TR>
</TABLE>
KEY TABLE

#---------End of Percentage Bases below 20 Table in 384 format---------------

#------Print Average Phred score in Qread-length Table in 384 format----------
print <<"ROW";
<BR>
<TABLE BORDER=0 COLS=2 WIDTH=600>
<TR>
<TD WIDTH=500>
<FONT SIZE=3><B><U>Average Phred Scores in <b>QRL</b>s for: <I>$Seqname</I></B></U></FONT><FONT SIZE=1> (Plate: $platename)</FONT>
<P>
ROW

print <<"AVG TABLE";
<TABLE BORDER=1 COLOR="000000" COLS=17 WIDTH=600>
AVG TABLE

$k = 5; $j = 81;
while ($k >= 0) {
$A = $start+$k; 

print <<"AVG TABLE";
<TR>
   <TD WIDTH=20><CENTER> $A </CENTER></TD>
AVG TABLE

   for($i=$j; $i <= $j+15; $i++) {
     $bcolor = $avgcolor[$i];
     $num = $avgvalue[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- bold & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
     $location = "trace1.cgi?filename=$Seqname&well=$i&start=$start";
     $jscript = "javascript:MessageScreen(\'$location\',\'$i\',500,300)\" style=\"text-decoration:none";
     print "<TD ALIGN=CENTER BGCOLOR=$bcolor WIDTH=20><A href=\"$jscript\"><FONT COLOR=$fcolor>$format</FONT></A></TD>\n";
   }
   $k--;
   $j -= 16;
   print "<TR>\n";
}
   
print <<"AVG TABLE";
<TR>   
   <TD WIDTH="20">&nbsp;</TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> A </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> B </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> C </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> D </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> E </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> F </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> G </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> H </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> I </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> J </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> K </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> L </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> M </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> N </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> O </FONT></CENTER></TD>
   <TD WIDTH="20"><CENTER><FONT COLOR=BLACK> P </FONT></CENTER></TD>
</TR>
</TABLE>
AVG TABLE

print <<"KEY TABLE";
</TD></TR>
<TR>
<TD ALIGN=LEFT WIDTH=400>
<TABLE BORDER=1 COLOR="000000" COLS=5 WIDTH=400>
<TR><TD ALIGN=CENTER><FONT SIZE=2> Color Key </FONT></TD>
    <TD BGCOLOR = #FF66FF ALIGN=CENTER><FONT SIZE=2> 0 - 20 </FONT></TD>
    <TD BGCOLOR = #33CCFF ALIGN=CENTER><FONT SIZE=2> 21 - 30 </FONT></TD>
    <TD BGCOLOR = #3333FF ALIGN=CENTER><FONT COLOR=WHITE SIZE=2> 31% - 40% </FONT></TD>
    <TD BGCOLOR = #0000BB ALIGN=CENTER><FONT COLOR=WHITE SIZE=2> 41+ </FONT></TD></TR>
</TABLE>
</TD>
</TR>
</TABLE>
KEY TABLE

#---------End of Average Phred score in Qread-length Table in 384 format---------------

#---------Print Capillary Value Table in 384 format-----------
if ($ABI == 1) {
$C = $start+5;
print <<"CAPILLARY";
<BR>
<FONT SIZE=3><B><U>Capillaries Used for: <I>$Seqname</I></B></U></FONT><FONT SIZE=1> (Plate: $platename)</FONT>
<P>
<TABLE BORDER=1 COLOR="000000" COLS=17 WIDTH=600>
CAPILLARY


$k = 5; $j = 81;
while ($k >= 0) {
$C = $start+$k; 
print <<"CAPILLARY";
<TR>
   <TD WIDTH=20><CENTER> $C </CENTER></TD>
CAPILLARY

   for($i=$j; $i <= $j+15; $i++) {
     $bcolor = $avgcolor[$i];
     $num = $Capillary[$i];
	
     $cflag = 0;
   
     if ($wC[$i] == 1) {  			# this is a control- underline & blink
        $format = "<BLINK><B> $num </B></BLINK>";
     } else {
     	$format = $num;
     }
     
     if(($bcolor eq "#3366FF") || ($bcolor eq "#3333FF") || ($bcolor eq "#0000BB")) {
	$fcolor = 'WHITE'; 
     }else {
        $fcolor = 'BLACK';       
     }
        
     print "<TD ALIGN=CENTER BGCOLOR=$bcolor WIDTH=20><FONT COLOR=$fcolor>$format</FONT></TD>\n";
   }
   $k--;
   $j -= 16;
   print "</TR>\n";
}

print <<"CAPILLARY";
<TR>   
   <TD WIDTH=20>&nbsp;</TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> A </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> B </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> C </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> D </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> E </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> F </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> G </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> H </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> I </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> J </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> K </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> L </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> M </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> N </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> O </FONT></CENTER></TD>
   <TD WIDTH=20><CENTER><FONT COLOR=BLACK> P </FONT></CENTER></TD>
</TR>
</TABLE>
CAPILLARY
}		# end of ABI Capillary Table
#-----End of Print Value Capillary Table in 384 format--------------------------------------
  }		# end else of 384 format display
}		# end of nodata == 0

#--------Footer for the HTML Page---------------
print <<"FOOTER PAGE";
<BR>
<HR>
<FONT SIZE=2>
For comments or questions regarding program scripts,  <BR>
please go to the <A HREF="http://wheat.pw.usda.gov/cgi-bin/mboard/feedback/list.cgi">Feedback Bulletin Board</A> 
<BR><HR>
Return to 
<A HREF="javascript:history.back()">Previous Page</A>
</FONT>

</BODY>
</HTML>
FOOTER PAGE
#--------End of Footer for HTML Page------------

#---------Subroutines---------------------------
sub defrow {
	$b = $_[0];
	$d = $_[1];

	if (($b eq "A")||($b eq "a")) {$d += 0;}	# A 1-12
	elsif(($b eq "B")||($b eq "b")) {$d += 12;}	# B 13-24
	elsif(($b eq "C")||($b eq "c")) {$d += 24;}	# C 25-36
	elsif(($b eq "D")||($b eq "d")) {$d += 36;}	# D 37-48
	elsif(($b eq "E")||($b eq "e")) {$d += 48;}	# E 49-60
	elsif(($b eq "F")||($b eq "f")) {$d += 60;}	# F 61-72
	elsif(($b eq "G")||($b eq "g")) {$d += 72;}	# G 73-84
	elsif(($b eq "H")||($b eq "h")) {$d += 84;}	# H 85-96
}

sub def384brow {
	$b = $_[0];			# well letter
	$d = $_[1];			# well number

	$offlet = &def384blet($b);
	$offnum = &def384bnum($d);
	$pos = $offlet + $offnum;	# well position
}

sub def384bnum {
	$d = $_[0];			# well number
	
	if (($d == 1) || ($d == 2)) { return 1; }
	elsif (($d == 3) || ($d == 4)) { return 2; }
	elsif (($d == 5) || ($d == 6)) { return 3; }
	elsif (($d == 7) || ($d == 8)) { return 4; }
	elsif (($d == 9) || ($d == 10)) { return 5; }	
	elsif (($d == 11) || ($d == 12)) { return 6; }
	elsif (($d == 13) || ($d == 14)) { return 7; }
	elsif (($d == 15) || ($d == 16)) { return 8; }
	elsif (($d == 17) || ($d == 18)) { return 9; }
	elsif (($d == 19) || ($d == 20)) { return 10; }
	elsif (($d == 21) || ($d == 22)) { return 11; }
	elsif (($d == 23) || ($d == 24)) { return 12; }
}

sub def384blet {
	$b = $_[0];			# well letter
	
	if (($b eq 'A') || ($b eq 'a')) { return 0; }
	elsif (($b eq 'B') || ($b eq 'b')) { return 0; }
	elsif (($b eq 'C') || ($b eq 'c')) { return 12; }
	elsif (($b eq 'D') || ($b eq 'd')) { return 12; }
	elsif (($b eq 'E') || ($b eq 'e')) { return 24; }
	elsif (($b eq 'F') || ($b eq 'f')) { return 24; }
	elsif (($b eq 'G') || ($b eq 'g')) { return 36; }
	elsif (($b eq 'H') || ($b eq 'h')) { return 36; }
	elsif (($b eq 'I') || ($b eq 'i')) { return 48; }
	elsif (($b eq 'J') || ($b eq 'j')) { return 48; }
	elsif (($b eq 'K') || ($b eq 'k')) { return 60; }
	elsif (($b eq 'L') || ($b eq 'l')) { return 60; }
	elsif (($b eq 'M') || ($b eq 'm')) { return 72; }
	elsif (($b eq 'N') || ($b eq 'n')) { return 72; }
	elsif (($b eq 'O') || ($b eq 'o')) { return 84; }
	elsif (($b eq 'P') || ($b eq 'p')) { return 84; }
}

sub def384row {
	$b = $_[0];			# well letter
	$d = $_[1];			# well number
	
	$offsetnum = $offsetlet = 0;	# initial the offset values
	
	$offsetnum = &def384num($d);	# assign offset according to num
	$offsetlet = &def384let($b);	# assign offset according to letter
	$pos = $offsetnum + $offsetlet; #  calculate the well index in 384-format
	
	return $pos;
}

sub def384num {
	$d = $_[0];			# well number
	
	if (($d == 1) || ($d == 7) || ($d == 13) || ($d == 19)) { return 0; }
	elsif (($d == 2) || ($d == 8) || ($d == 14) || ($d == 20)) { return 16; }
	elsif (($d == 3) || ($d == 9) || ($d == 15) || ($d == 21)) { return 32; }
	elsif (($d == 4) || ($d == 10) || ($d == 16) || ($d == 22)) { return 48; }
	elsif (($d == 5) || ($d == 11) || ($d == 17) || ($d == 23)) { return 64; }	
	elsif (($d == 6) || ($d == 12) || ($d == 18) || ($d == 24)) { return 80; }
}

sub def384let {
	$b = $_[0];			# well letter
	
	if (($b eq 'A') || ($b eq 'a')) { return 1; }
	elsif (($b eq 'B') || ($b eq 'b')) { return 2; }
	elsif (($b eq 'C') || ($b eq 'c')) { return 3; }
	elsif (($b eq 'D') || ($b eq 'd')) { return 4; }
	elsif (($b eq 'E') || ($b eq 'e')) { return 5; }
	elsif (($b eq 'F') || ($b eq 'f')) { return 6; }
	elsif (($b eq 'G') || ($b eq 'g')) { return 7; }
	elsif (($b eq 'H') || ($b eq 'h')) { return 8; }
	elsif (($b eq 'I') || ($b eq 'i')) { return 9; }
	elsif (($b eq 'J') || ($b eq 'j')) { return 10; }
	elsif (($b eq 'K') || ($b eq 'k')) { return 11; }
	elsif (($b eq 'L') || ($b eq 'l')) { return 12; }
	elsif (($b eq 'M') || ($b eq 'm')) { return 13; }
	elsif (($b eq 'N') || ($b eq 'n')) { return 14; }
	elsif (($b eq 'O') || ($b eq 'o')) { return 15; }
	elsif (($b eq 'P') || ($b eq 'p')) { return 16; }
}

sub defcolor {
	$a = $_[0];

	if(($a >= 0.0) && ($a <= 0.1))  # Define Color
                {$value = "#FF00FF";}		
	elsif (($a > 0.1) && ($a <= 0.2))
		{$value = "#FF66FF";}		
	elsif (($a > 0.2) && ($a <= 0.3))
		{$value = "#FF99FF";}
	elsif (($a > 0.3) && ($a <= 0.4))
		{$value = "#FFCCFF";}
	elsif (($a > 0.4) && ($a <= 0.5))
		{$value = "#66CCFF";}
	elsif (($a > 0.5) && ($a <= 0.6))
		{$value = "#33CCFF";}
	elsif (($a > 0.6) && ($a <= 0.7))
		{$value = "#3399FF";}
	elsif (($a > 0.7) && ($a <= 0.8))
		{$value = "#3366FF";}
	elsif (($a > 0.8) && ($a <= 0.9))
		{$value = "#3333FF";}
	elsif (($a > 0.9) && ($a <= 1.0))
		{$value = "#0000BB";}
}

sub deGcolour {
	$c = $_[0];

	if(($c >= 0) && ($c <= 150))		# Define Good Value color
		{$goodcolour = "#FF66FF";}
	elsif (($c >= 151) && ($c <= 300))
		{$goodcolour = "#FF99FF";}
	elsif (($c >= 301) && ($c <= 450))
		{$goodcolour = "#33CCFF";}
	elsif (($c >= 451) && ($c <= 600))
		{$goodcolour = "#3333FF";}
	elsif ($c >= 601)
		{$goodcolour = "#0000BB";}
}

sub defbadcolor {
	$b = $_[0];

	if(($b >= 0) && ($b <= 20))		# Define Bad Value color
		{$badcolor = "#0000BB";}
	elsif (($b >= 21) && ($b <= 30))
		{$badcolor = "#33CCFF";}
	elsif (($b >= 31) && ($b <= 40))
		{$badcolor = "#FF99FF";}
	elsif ($b >=41)
		{$badcolor = "#FF66FF";}
}

sub defavgcolor {
	$c = $_[0];

	if(($c >= 0) && ($c <= 20))		# Define Average Value color
		{$avgcolor = "#FF66FF";}
	elsif (($c >= 21) && ($c <= 30))
		{$avgcolor = "#33CCFF";}
	elsif (($c >= 31) && ($c <= 40))
		{$avgcolor = "#3333FF";}
	elsif ($c >=41)
		{$avgcolor = "#0000BB";}
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
