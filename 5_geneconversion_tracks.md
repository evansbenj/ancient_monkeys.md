# Identifying gene conversion tracks.

I wrote a script to do this for reciprocally monophyletic clades (gets_geneconversion_tracks__.pl):

```
#!/usr/bin/perl 
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use Array::Utils qw(:all);


# TO execute this program type this:
# ./gets_geneconversion_tracks.pl arg1 arg2
# where arg1 and arg2 are input and output file names, respectively
# This program reads a file with an outgroup sequences and sequences from two species
# the first line is the number of bp, and the number of sequences from the first population 
# and then the second
# the program then identifies sites that are not consistent with monophyly of each population
# for each of these sites it will figure out which samples match each other
# and it will look for runs of sites where the same individuals are similarly discordant.

# the first sequence should be an outgroup

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

#### Prepare the input file with sequences

unless (open DATAINPUT, $inputfile) {
	print "Can not find the data input file!\n";
	exit;
}

#my $outputfile = "geneconversion_tracks.out";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";



###################
my $linenumber=0;
my %hash;
my @numbers;
my @temp;
my @temp1;
my $m;
my $n;
my $r;
my $z;
my $q;
my $invariant=0;
my $invariantpop1=0;
my $invariantpop2=0;
my $ancestralpop1;
my $ancestralpop2;
my @column;
my @column_pop1;
my @column_pop2;
my @uniq_chars;
my @uniq_chars_pop1;
my @uniq_chars_pop2;
my %varhash;
my $counter_pop1=0;
my $counter_pop2=0;
my %group1; # this is a hash with bp as position that has the individuals in group1
my %group2;
my %snake1; # this is a hash with bp as position that has the individuals in snake1
my %snake2;


while ( my $line = <DATAINPUT>) {
## set first line equal to gene name

	if ($linenumber == 0)  { # this is the first line
		@numbers = split(/\s+/,$line);
		print  "The number of bases is: $numbers[0] \n";
		print  "The number of seqs in pop1 is $numbers[1] \n";
		print  "The number of seqs in pop2 is $numbers[2] \n";
				$linenumber += 1;
	}
	elsif (($linenumber ne 0) && ($linenumber <= ($numbers[1]+$numbers[2]+1))) { 
	# this is the line containing the sequences (the first is the outgroup)
		@temp = split(/\s+/,$line);
		$hash{$linenumber-1}[0]=$temp[0];
		@temp1=split('',$temp[1]);
		for ($m=0; $m <= $#temp1; $m++){
			$hash{$linenumber-1}[$m+1]=$temp1[$m];
		}
		@temp=();  #this should erase the contents of the array
		@temp1=();  #this should erase the contents of the array
				$linenumber += 1;
	} 
}	#end while

# now we have a hash where with a name in [0] and bp in all of the rest.

#cycle through each position and categorize it
for ($n=1; $n <= $numbers[0]; $n++){
		#first test if the site is invariant in ingroup
		for ($r=1; $r <= ($numbers[1]+$numbers[2]); $r++){
			push(@column, $hash{$r}[$n]);
		} #end $r
		#print "@column\n";
		@uniq_chars = uniq @column;
		if ($#uniq_chars == 1){ # we are going to ignore positions with more than 2 variants (brings it down to 706 from 829)
			# this is a variable site
			# now check if it is one that violates the monophyly of each pop
			#first test if the site is invariant in ingroup
			@column_pop1 = @column[0..($numbers[1]-1)];
			@column_pop2 = @column[$numbers[1]..($numbers[1]+$numbers[2]-1)];
			@uniq_chars_pop1 = uniq @column_pop1;
			@uniq_chars_pop2 = uniq @column_pop2;
			if(($#uniq_chars_pop1==1)&&($#uniq_chars_pop2==1)&& (sort (@uniq_chars_pop1) ~~ sort (@uniq_chars_pop2))){ 
				# both groups have the same SNPs
				# now count to see if both groups have at least 2 individuals with each SNP
				$counter_pop1=0;
				$counter_pop2=0;
				for ($q=0; $q <= $#column_pop1; $q++){
					if($column_pop1[$q] eq $uniq_chars_pop1[0]){
						$counter_pop1+=1
					}
				}	
				for ($q=0; $q <= $#column_pop2; $q++){
					if($column_pop2[$q] eq $uniq_chars_pop2[0]){
						$counter_pop2+=1
					}
				}	
				if(($counter_pop1 != 1)&&($counter_pop1 != $#column_pop1)&&($counter_pop2 != 1)&&($counter_pop2 != $#column_pop2)){
					# each group does not have a singleton SNP
					print "$n @column_pop1 ZZZ @column_pop2 \n";
					# add these sites to a new data array; this will be used to identify tracks
					for ($z=0; $z <= ($numbers[1]+$numbers[2]); $z++){
						$varhash{$n}[$z] = $column[$z];
					}
				}	
			}	
		}
		@column=();
}

# ok now we have %varhash loaded with the interesting positions
# now let's run through this and look for gene conversion tracks.
# these are defined as (1) consecutive variable positions that (2) have at least two taxon from each clade 
# that violates monophyly and where (3) consecutive sites have at least two taxon from both clades and 
# (4) tracks are not interrupted by other sites that meet criterion 2

my $keycounter=0;
my $key;
my @positionz;
my $counter1=0;
my $counter2=0;
# make two arrays with the number of individuals from each clade
my @gp1_individuals_array;
my @gp2_individuals_array;

@gp1_individuals_array = (0..($numbers[1]-1));
@gp2_individuals_array = ($numbers[1]..($numbers[1]+$numbers[2]-1));

# cycle through the variable positions and load group hashes with groups for each interesting position
foreach my $key (sort { $a <=> $b } keys(%varhash) ){
	# for this column in the alignment, define the two groups using the first position
	$counter1=0;
	$counter2=0;
	$group1{$key}[$counter1]=0;
	$counter1+=1;
	for ($z=1; $z < ($numbers[1]+$numbers[2]); $z++){
		if($varhash{$key}[0] eq $varhash{$key}[$z]){
			$group1{$key}[$counter1]=$z;
			$counter1+=1;
		}
		else{
			$group2{$key}[$counter2]=$z;
			$counter2+=1;
		}	
	}	
}
# now we have two hashes of arrays that each contain the individuals that match 
# for each position of the alignment (which is the key of the hashes)
$counter1=0;
$counter2=0;
my $previous_key=0;
my @isect1=();
my @isect2=();
my @isect3=();
my @isect4=();
my @isect5=();
my @isect6=();
my @isect7=();
my @isect8=();
my @isect9=();
my @isect10=();
my @isect11=();
my @isect12=();

# now cycle through group hash and calculate intersections between adjacent positions 
# and assign snake hash based on intersections and shared SNPs in each clade

foreach my $key (sort { $a <=> $b } keys(%group1) ){
	if($previous_key != 0){
		@isect1=();
		@isect2=();
		@isect3=();
		@isect4=();
		@isect5=();
		@isect6=();
		@isect7=();
		@isect8=();
		@isect9=();
		@isect10=();
		@isect11=();
		@isect12=();
		#print "previous_gp_hey $previous_key @{ $group1{$previous_key} } GGG @{ $group2{$previous_key} } \n";
		#print "gp_hey $key @{ $group1{$key} } GGG @{ $group2{$key} } \n";
		#print "hello $key @{ $group1{$key} } YYYYYYYYY @{ $group2{$key} }\n";
		@isect1 = intersect(@{ $group1{$key} }, @{ $group1{$previous_key} });
		@isect2 = intersect(@{ $group2{$key} }, @{ $group1{$previous_key} });
		@isect3 = intersect(@{ $group1{$key} }, @{ $group2{$previous_key} });
		@isect4 = intersect(@{ $group2{$key} }, @{ $group2{$previous_key} });
		#print "hello $key @isect1 YYYYYYYYY @isect2 ZZZZZ @isect3 AAAA @isect4\n";
		# now select the two intersections that are consistent with gene conversion and add the larger one to the snake
		# if there are two.  In order for them to be consistent with gene conversion, both intersections have to have
		# two members from each of both clades
		@isect5 = intersect(@isect1, @gp1_individuals_array);
		@isect6 = intersect(@isect2, @gp1_individuals_array);
		@isect7 = intersect(@isect3, @gp1_individuals_array);
		@isect8 = intersect(@isect4, @gp1_individuals_array);
		@isect9 = intersect(@isect1, @gp2_individuals_array);
		@isect10 = intersect(@isect2, @gp2_individuals_array);
		@isect11 = intersect(@isect3, @gp2_individuals_array);
		@isect12 = intersect(@isect4, @gp2_individuals_array);
		print "goodby $key @isect5 YYYYYYYYY @isect6 ZZZZZ @isect7 AAAA @isect8 FFF @isect9 YYYYYYYYY @isect10 ZZZZZ @isect11 AAAA @isect12\n";
		# now check if both tracks have two members from each group
		print $key," DDD ",$#isect5," DDD ",$#isect9," DDD ",$#isect8," DDD ",$#isect12,"\n";
		if(($#isect5 >= 1)&&($#isect9 >= 1)&&($#isect8 >= 1)&&($#isect12 >= 1)){
			push(@{$snake1{$key}}, @isect1);
			push(@{$snake2{$key}}, @isect4);
			if(exists $snake1{$previous_key-0.5} ){
				push(@{$snake1{$previous_key}}, @isect1); # this puts a subset of the (appropriate) group in the beginning of the first position of
														  # the track based on how it matches the next position
				push(@{$snake2{$previous_key}}, @isect4);
			}
			#print "wacko1 ",@{ $snake1{$key} },"\n";
		}
		elsif(($#isect7 > 1)&&($#isect11 > 1)&&($#isect6 > 1)&&($#isect10 > 1)){
			push(@{$snake1{$key}}, @isect3);
			push(@{$snake2{$key}}, @isect2);
			if(exists $snake1{$previous_key-0.5} ){
				push(@{$snake1{$previous_key}}, @isect3);
				push(@{$snake2{$previous_key}}, @isect2);
			}
			#print "wacko2 ",@{ $snake1{$key} },"\n";
		}
		else{
			push(@{$snake1{$key-0.5}}, "no_extension");
			push(@{$snake2{$key-0.5}}, "no_extension");
			#print "wacko3 ",@{ $snake1{$key} },"\n";
		}
		$previous_key=$key;
	}
	else{
		$previous_key=$key;
		push(@{$snake1{$key}}, "first_position");
		push(@{$snake2{$key}}, "first_position");
		#print "wacko4 ",@{ $snake1{$key} },"\n";
	}	
}

foreach $key (sort { $a <=> $b } keys(%snake1) ){
  	print OUTFILE "$key @{ $snake1{$key} } \|\|\| @{ $snake2{$key} }\n";
  	print "$key @{ $snake1{$key} } \|\|\| @{ $snake2{$key} }\n";
}	

# now print to screen the ranges of geneconv tracks; assume a minimum length of $minlength positions
my $minlength=3;
$counter1=0;
my $begin;

foreach $key (sort { $a <=> $b } keys(%snake1) ){
	#print printf("%.0f",$key)," ",$key,"\n";
	my $key2=sprintf("%.0f",$key);
	if($key2 != $key){
		# this is the end of a track
		if($counter1>($minlength)){
			print $begin," - ",$previous_key," ";
		}
		$counter1=0;
	}
	elsif($counter1==1){
		$begin=$key;
	}
	elsif($snake1{$key}[0]  =~ 'first'){
		# this is the beginning of a track
		$begin=$key;
	}
	$previous_key=$key;
	$counter1+=1;
}
if($counter1>($minlength)){
	print $begin," - ",$previous_key," ";
}
print "\n";



```


I also wrote a perl script to parse the output of RDP (interprets_RDP_output.pl):
```
#!/usr/bin/perl 
use strict;
use warnings;
no warnings qw(uninitialized);
use List::MoreUtils 'true';
use Array::Utils qw(:all);



# TO execute this program type this:
# ./interprets_RDP_output.pl arg1 arg2
# where arg1 and arg2 are input and output file names, respectively
# This script reads a file from the program RDP (recombination detection program)
# the script then will sort the output into ranges for each method based on criteria I will set.
# for now I will get start and correspinding end points of ranges
# that involve nemestrina and sulawesi
# to do the same for the fascicularis output, you need to comment out the test for nem and sula
# and uncomment the test for fasc clades A and B

# the first sequence should be an outgroup

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

#### Prepare the input file with sequences

unless (open DATAINPUT, $inputfile) {
	print "Can not find the data input file!\n";
	exit;
}

#my $outputfile = "geneconversion_tracks.out";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";



###################
my $linenumber=0;
my @temp;
my $m;
my @recomb_array;
my @major_parent;
my @minor_parent;
my $start;
my $end;
my $nem_found1;
my $sula_found1;
my $nem_found2;
my $sula_found2;
my $nem_found3;
my $sula_found3;
my $previous_temp_0="1";
my $RDPflag;
my $geneconvflag;
my $bootscanflag;
my $maxchiflag;
my $chimaeraflag;
my $sisscanflag;
my $phyloproflag;
my $lardflag;
my $threeseqflag;
my @RDP_start=();
my @geneconv_start=();
my @bootscan_start=();
my @maxchi_start=();
my @chimaera_start=();
my @sisscan_start=();
my @phylopro_start=();
my @lard_start=();
my @threeseq_start=();
my @RDP_stop=();
my @geneconv_stop=();
my @bootscan_stop=();
my @maxchi_stop=();
my @chimaera_stop=();
my @sisscan_stop=();
my @phylopro_stop=();
my @lard_stop=();
my @threeseq_stop=();


while ( my $line = <DATAINPUT>) {
	if (($linenumber == 0)||($linenumber == 1)||($linenumber == 2))  { # these are only comments
		$linenumber += 1;
	}
	else { # read in the data
		chomp $line; # get rid of carriage return
		if(length($line)>0){ # this is supposed to avoid dealing with blank lines
			@temp = split(',',$line);
			$temp[0]=~ s/^\s+|\s+$//g; # get rid of spaces
			if(@temp){ # this is supposed to deal with the blank lines
				s/\*// for @temp; # get rid of asterisks
				if($temp[0] eq $previous_temp_0){
					if($temp[0] eq "1"){ # this records the start and stop of the first line
											# for the other lines this is recorded later
						$start=$temp[2];
						$end=$temp[3];					
					}
					# this is the same region so add names to the current track
					push(@recomb_array,$temp[8]);
					push(@major_parent,$temp[9]);
					push(@minor_parent,$temp[10]);
				}
				else{
					# this is a new region so process the previous track
					if((defined $temp[2])&&(defined $temp[3])&&(@recomb_array)&&(@major_parent)&&(@minor_parent)){ # check if the start and stop are both numbers
						# this is the test for nem/Sula
						$nem_found1 = true { /nem/ } @recomb_array;
						$sula_found1 = true { /heck|tonk/ } @recomb_array;
						$nem_found2 = true { /nem/ } @major_parent;
						$sula_found2 = true { /heck|tonk/ } @major_parent;
						$nem_found3 = true { /nem/ } @minor_parent;
						$sula_found3 = true { /heck|tonk|maur|nig|bru|tog/ } @minor_parent;
						# this is the end of the test for nem/sula

						# this is the test for fasc A/B
						#$nem_found1 = true { /A_/ } @recomb_array;
						#$sula_found1 = true { /B_/ } @recomb_array;
						#$nem_found2 = true { /A_/ } @major_parent;
						#$sula_found2 = true { /B_/ } @major_parent;
						#$nem_found3 = true { /A_/ } @minor_parent;
						#$sula_found3 = true { /B_/ } @minor_parent;
						# this is the end of the test for fasc A/B
						print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
						if(
							(($nem_found1>=1)&&($sula_found1>=1))||(($nem_found1>=1)&&($sula_found2>=1))||(($nem_found1>=1)&&($sula_found3>=1)) ||
							(($nem_found2>=1)&&($sula_found1>=1))||(($nem_found2>=1)&&($sula_found2>=1))||(($nem_found2>=1)&&($sula_found3>=1)) ||
							(($nem_found3>=1)&&($sula_found1>=1))||(($nem_found3>=1)&&($sula_found2>=1))||(($nem_found3>=1)&&($sula_found3>=1))
							){
								if(($RDPflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@RDP_start,$start);
									push(@RDP_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($geneconvflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@geneconv_start,$start);
									push(@geneconv_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($bootscanflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@geneconv_start,$start);
									push(@geneconv_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($geneconvflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@bootscan_start,$start);
									push(@bootscan_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($maxchiflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@maxchi_start,$start);
									push(@maxchi_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($chimaeraflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@chimaera_start,$start);
									push(@chimaera_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($sisscanflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@sisscan_start,$start);
									push(@sisscan_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($phyloproflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@phylopro_start,$start);
									push(@phylopro_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($lardflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@lard_start,$start);
									push(@lard_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
								if(($threeseqflag ne 'NS')&&($start ne "")&&($end ne "")){
									push(@threeseq_start,$start);
									push(@threeseq_stop,$end);	
									#print "$nem_found1 $sula_found1 $linenumber you $previous_temp_0 GGG $temp[0] YYY @major_parent and @minor_parent start $start end $end\n";
								}
						}
						if(($temp[2] ne "")&&($temp[3] ne "")){
							$start=$temp[2];
							$end=$temp[3];
						}
						else{
							print "yo ",$temp[2]," ",$temp[2],"\n";
						}
						# and now delete arrays and start a new one
						$previous_temp_0=$temp[0];
						@recomb_array=();
						@major_parent=();
						@minor_parent=();
						push(@recomb_array,$temp[8]);
						push(@major_parent,$temp[9]);
						push(@minor_parent,$temp[10]);
						$RDPflag=$temp[11];
						$geneconvflag=$temp[12];
						$bootscanflag=$temp[13];
						$maxchiflag=$temp[14];
						$chimaeraflag=$temp[15];
						$sisscanflag=$temp[16];
						$phyloproflag=$temp[17];
						$lardflag=$temp[18];
						$threeseqflag=$temp[19];
					}
				}
			$linenumber += 1;
			}
		}	
	}	
	# now process the last array
}	#end while

print "RDP start ",$#RDP_start+1," ranges\n";
print OUTFILE "RDP start ",$#RDP_start+1," ranges\n";
foreach(@RDP_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "RDP end\n";
print OUTFILE "RDP end\n";
foreach(@RDP_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "geneconv start ",$#geneconv_start+1," ranges\n";
print OUTFILE "geneconv start ",$#geneconv_start+1," ranges\n";
foreach(@geneconv_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "geneconv end\n";
print OUTFILE "geneconv end\n";
foreach(@geneconv_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "bootscan start ",$#bootscan_start+1," ranges\n";
print OUTFILE "bootscan start ",$#bootscan_start+1," ranges\n";
foreach(@bootscan_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "bootscan end\n";
print OUTFILE "bootscan end\n";
foreach(@bootscan_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "maxchi start ", $#maxchi_start+1," ranges\n";
print OUTFILE "maxchi start ", $#maxchi_start+1," ranges\n";
foreach(@maxchi_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "maxchi end\n";
print OUTFILE "maxchi end\n";
foreach(@maxchi_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "chimaera start ",$#chimaera_start+1, " ranges\n";
print OUTFILE "chimaera start ",$#chimaera_start+1, " ranges\n";
foreach(@chimaera_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "chimaera end\n";
print OUTFILE "chimaera end\n";
foreach(@chimaera_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "sisscan start ",$#sisscan_start+1 , " ranges\n";
print OUTFILE "sisscan start ",$#sisscan_start+1 , " ranges\n";
foreach(@sisscan_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "sisscan end\n";
print OUTFILE "sisscan end\n";
foreach(@sisscan_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "phylopro start ",$#phylopro_start+1," ranges\n";
print OUTFILE "phylopro start ",$#phylopro_start+1," ranges\n";
foreach(@phylopro_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "phylopro end\n";
print OUTFILE "phylopro end\n";
foreach(@phylopro_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "lard start ",$#lard_start+1, " ranges\n";
print OUTFILE "lard start ",$#lard_start+1, " ranges\n";
foreach(@lard_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "lard end\n";
print OUTFILE "lard end\n";
foreach(@lard_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "threeseq start ",$#threeseq_start+1, " ranges\n";
print OUTFILE "threeseq start ",$#threeseq_start+1, " ranges\n";
foreach(@threeseq_start){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";
print "threeseq end\n";
print OUTFILE "threeseq end\n";
foreach(@threeseq_stop){
	print $_,",";
	print OUTFILE $_,",";
}
print "\n";
print OUTFILE "\n";

```

And I wrote an R script to plot the results (range_plot.R):
```
# This R script will (hopefully) allow me to present gene conversion tracks
# from several programs simultaneously

# Install GenomicRanges package if this is the first time
#if (!require("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")

# Install GenomicRanges package if this is the first time
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

# this is an old verision
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap", update = TRUE, ask = FALSE)


# install for first time
#install.packages("UpSetR")
#devtools::install_github("hms-dbmi/UpSetR")
library(circlize)
library(GenomicRanges)
library(UpSetR)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)

# 
# information/ideas from here: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

# make GRange variables for each method
# First let's add mine (reciprocal_snakes)

reciprocal_snakes <- GRanges(
  seqnames = Rle(c("snak"), c(6)),
  ranges = IRanges(start = c(5390,6058,9110,10015,13364,16456), 
                   end = c(5906,9088,9622,11024,14104,16612)
  )
)
#reciprocal_snakes

rdp <- GRanges(
  seqnames = Rle(c("rdp"), c(34)),
  ranges = IRanges(start = c(8154, 7684, 7463, 7864, 7670, 7214, 6684, 7694, 8272, 16456, 6606, 10174, 6531, 10015, 15770, 7544, 5744, 7008, 6906, 8197, 11775, 5546, 7477, 9915, 8272, 10189, 15749, 6796, 5554, 7298, 5510, 5390, 3912, 10394), 
                   end = c(8606, 7960, 8153, 8196, 8017, 7470, 7196, 7876, 8554, 16612, 7472, 10493, 8832, 10476, 15872, 8141, 5918, 7425, 7425, 8288, 12083, 5906, 7611, 10588, 8832, 11024, 15875, 7213, 5830, 7552, 5743, 5556, 4186, 10588)
  )
)
#rdp

geneconv <- GRanges(
  seqnames = Rle(c("genc"), c(79)),
  ranges = IRanges(start = c(8154, 8154, 7684, 7684, 7463, 7463, 7864, 7864, 7670, 7670, 7214, 7214, 6684, 6684, 7694, 7694, 8272, 8272, 16456, 16456, 6606, 6606, 10174, 10174, 6531, 6531, 10015, 10015, 15770, 15770, 7544, 7544, 5744, 5744, 7008, 7008, 6906, 6906, 8197, 8197, 11775, 11775, 5546, 5546, 7477, 7477, 9915, 9915, 8272, 8272, 8272, 10189, 10189, 5744, 15749, 15749, 649, 649, 6796, 6796, 5554, 5554, 5661, 5661, 7426, 11408, 11408, 7298, 7298, 9977, 9977, 5510, 5510, 8273, 5390, 5390, 3912, 3912, 10394),
                   end = c(8606, 8606, 7960, 7960, 8153, 8153, 8196, 8196, 8017, 8017, 7470, 7470, 7196, 7196, 7876, 7876, 8554, 8554, 16612, 16612, 7472, 7472, 10493, 10493, 8832, 8832, 10476, 10476, 15872, 15872, 8141, 8141, 5918, 5918, 7425, 7425, 7425, 7425, 8288, 8288, 12083, 12083, 5906, 5906, 7611, 7611, 10588, 10588, 8832, 8832, 8488, 11024, 11024, 5906, 15875, 15875, 757, 757, 7213, 7213, 5830, 5830, 5920, 5920, 8271, 11616, 11616, 7552, 7552, 10042, 10042, 5743, 5743, 8393, 5556, 5556, 4186, 4186, 10588)
  )
)
#geneconv

bootscan <- GRanges(
  seqnames = Rle(c("boot"), c(41)),
  ranges = IRanges(start = c(8154, 7684, 7463, 7864, 7670, 7214, 6684, 7694, 8272, 16456, 6606, 10174, 6531, 10015, 15770, 7544, 5744, 7008, 6906, 8197, 11775, 5546, 7477, 9915, 8272, 8272, 10189, 5744, 15749, 649, 6796, 5554, 5661, 7426, 11408, 7298, 9977, 5510, 8273, 5390, 3912),
                   end = c(8606, 7960, 8153, 8196, 8017, 7470, 7196, 7876, 8554, 16612, 7472, 10493, 8832, 10476, 15872, 8141, 5918, 7425, 7425, 8288, 12083, 5906, 7611, 10588, 8832, 8488, 11024, 5906, 15875, 757, 7213, 5830, 5920, 8271, 11616, 7552, 10042, 5743, 8393, 5556, 4186)
  )
)
#bootscan
maxchi <- GRanges(
  seqnames = Rle(c("maxX"), c(33)),
  ranges = IRanges(start = c(8154, 7684, 7463, 7864, 7670, 7214, 6684, 7694, 8272, 16456, 6606, 10174, 6531, 10015, 15770, 7544, 5744, 7008, 6906, 8197, 11775, 5546, 7477, 9915, 8272, 10189, 15749, 6796, 5554, 7426, 7298, 8273, 5390),
            end = c(8606, 7960, 8153, 8196, 8017, 7470, 7196, 7876, 8554, 16612, 7472, 10493, 8832, 10476, 15872, 8141, 5918, 7425, 7425, 8288, 12083, 5906, 7611, 10588, 8832, 11024, 15875, 7213, 5830, 8271, 7552, 8393, 5556)
  )
)

chimaera <- GRanges(
  seqnames = Rle(c("chim"), c(31)),
  ranges = IRanges(start = c(8154, 7684, 7463, 7864, 7670, 7214, 6684, 7694, 8272, 16456, 6606, 10174, 6531, 10015, 15770, 7544, 5744, 7008, 6906, 8197, 11775, 5546, 7477, 9915, 8272, 10189, 6796, 7426, 7298, 8273, 5390),
                   end = c(8606, 7960, 8153, 8196, 8017, 7470, 7196, 7876, 8554, 16612, 7472, 10493, 8832, 10476, 15872, 8141, 5918, 7425, 7425, 8288, 12083, 5906, 7611, 10588, 8832, 11024, 7213, 8271, 7552, 8393, 5556)
  )
)
sisscan <- GRanges(
  seqnames = Rle(c("siss"), c(23)),
  ranges = IRanges(start = c(8154, 7684, 7463, 7864, 7670, 6684, 16456, 6606, 10174, 6531, 10015, 7008, 6906, 8197, 5546, 7477, 8272, 10189, 15749, 7426, 11408, 5510, 5390),
                   end = c(8606, 7960, 8153, 8196, 8017, 7196, 16612, 7472, 10493, 8832, 10476, 7425, 7425, 8288, 5906, 7611, 8832, 11024, 15875, 8271, 11616, 5743, 5556)
  )
)

threeseq <- GRanges(
  seqnames = Rle(c("3seq"), c(31)),
  ranges = IRanges(start = c(8154, 7684, 7463, 7864, 7670, 7214, 6684, 7694, 8272, 16456, 6606, 10174, 6531, 10015, 15770, 5744, 7008, 6906, 8197, 11775, 7477, 9915, 8272, 8272, 5744, 6796, 5554, 5661, 7426, 5390, 3912),
                   end = c(8606, 7960, 8153, 8196, 8017, 7470, 7196, 7876, 8554, 16612, 7472, 10493, 8832, 10476, 15872, 5918, 7425, 7425, 8288, 12083, 7611, 10588, 8832, 8488, 5906, 7213, 5830, 5920, 8271, 5556, 4186)
  )
)

phylopro <- GRanges(
  seqnames = Rle(c("phylopro"), c(0)),
  ranges = IRanges(start = c(),
                   end = c()
  )
)
lard <- GRanges(
  seqnames = Rle(c("lard"), c(0)),
  ranges = IRanges(start = c(),
                   end = c()
  )
)

df_reciprocal_snakes <- as.data.frame(reciprocal_snakes)
df_rdp <- as.data.frame(rdp)
df_geneconv <- as.data.frame(geneconv)
df_bootscan <- as.data.frame(bootscan)
df_maxchi <- as.data.frame(maxchi)
df_chimaera <- as.data.frame(chimaera)
df_sisscan <- as.data.frame(sisscan)
df_threeseq <- as.data.frame(threeseq)
df_phylopro <- as.data.frame(phylopro)
df_lard <- as.data.frame(lard)
all_df<-rbind(df_reciprocal_snakes,df_rdp,df_geneconv,df_bootscan,df_maxchi,df_chimaera,df_sisscan,df_threeseq,df_phylopro,df_lard)

ggplot(all_df, aes(x=seqnames))+
  geom_linerange(aes(ymin=start,ymax=end),linetype=1,color="red")+
  geom_point(aes(y=start),size=0.25,color="red")+
  geom_point(aes(y=end),size=0.25,color="red")+
  theme_bw() 



### Below not used

#lt = list(reciprocal_snakes,rdp,geneconv)

#m3 = make_comb_mat(lt, mode = "union")

#getNamespaceExports("ComplexHeatmap")


```
