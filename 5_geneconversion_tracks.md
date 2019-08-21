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
		#print "goodby $key @isect5 YYYYYYYYY @isect6 ZZZZZ @isect7 AAAA @isect8 FFF @isect9 YYYYYYYYY @isect10 ZZZZZ @isect11 AAAA @isect12\n";
		# now check if both tracks have two members from each group
		#print $key," DDD ",$#isect5," DDD ",$#isect9," DDD ",$#isect8," DDD ",$#isect12,"\n";
		if(($#isect5 >= 1)&&($#isect9 >= 1)&&($#isect8 >= 1)&&($#isect12 >= 1)){
			push(@{$snake1{$key}}, @isect1);
			push(@{$snake2{$key}}, @isect4);
			#print "wacko1 ",@{ $snake1{$key} },"\n";
		}
		elsif(($#isect7 > 1)&&($#isect11 > 1)&&($#isect6 > 1)&&($#isect10 > 1)){
			push(@{$snake1{$key}}, @isect3);
			push(@{$snake2{$key}}, @isect2);
			#print "wacko2 ",@{ $snake1{$key} },"\n";
		}
		else{
			push(@{$snake1{$key}}, "no_extension");
			push(@{$snake2{$key}}, "no_extension");
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

```
