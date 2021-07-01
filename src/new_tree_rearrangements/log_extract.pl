use strict;
use warnings;
my @logs=grep { !/.*pb$/ } glob("*");
for my $name ( @logs){
    $name=~/(\d+)\.(.*)/;
    print "$1\t$2\t";
    open(my $fh,"<",$name);
    my $is_first=1;
    while(<$fh>){
        if($_=~/after state reassignment:(\d+)/){
            print "$1\t";
        }
        if($_=~/Will search ([\d.]+) of nodes/ && $is_first){
            print "$1\t";
            $is_first=0;
        }
        if($_=~/Final Parsimony score ([\d.]+)/){
            if($is_first){
                print "1\t";
            }
            print "$1\n";
        }
    }
}