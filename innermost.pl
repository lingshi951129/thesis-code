#!/usr/bin/perl
$prefix=$ARGV[0];
$in=$prefix.'multi_domain.txt';
$out1=$prefix.'multi_inner.txt';
$out2=$prefix.'delete.txt';
open GM,"$in";
$tab="\t";
$line=<GM>;
while($line=<GM>){
    if(!($line=~/^\*/)){
        chomp $line;
        @info=split(/\t/,$line);
        $gene=$info[0];
        if(exists $head{$gene} && $mark{$gene} !=1){ 
            if($info[4]>=$head{$gene} && $info[5]<=$tail{$gene}){
                $head{$gene}=$info[4];
                $tail{$gene}=$info[5];
                $tad{$gene}=$info[3];
            }
            elsif($info[4]<=$head{$gene} && $info[5]>=$tail{$gene}){
                next;
            }
            else{
                $mark{$gene}=1;
                delete $head{$gene};
                delete $tail{$gene};
                delete $tad{$gene};
                delete $chr{$gene};
                delete $pos{$gene};
                $delete{$gene}=$info[0].$tab.$info[1].$tab.$info[2];
            }
        }
        elsif(!(exists $delete{$gene}) && $mark{$gene} !=1){
            $head{$gene}=$info[4];
            $tail{$gene}=$info[5];
            $chr{$gene}=$info[1];
            $pos{$gene}=$info[2];
            $tad{$gene}=$info[3];
        }
    }
}
open GMOUT,">$out1";
print GMOUT "gene\tchr\tpos\ttadmark\ttadhead\ttadtail\n";
foreach my $key(keys %head){
    if($mark{$key} !=1){
        print GMOUT "$key\t$chr{$key}\t$pos{$key}\t$tad{$key}\t$head{$key}\t$tail{$key}\n";
    }
}
open GMDEL,">$out2";
print GMDEL "gene\tchr\tpos\n";
foreach my $key(keys %delete){
    print GMDEL "$delete{$key}\n";
}
