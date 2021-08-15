#!/usr/bin/perl
use List::Util qw(min max);
open SAME,"./same_compare_1000kb_TPM_flt.txt";
open DIF,"./dif_compare_1000kb_TPM_flt.txt";
$line=<SAME>;
while($line=<SAME>){
    chomp $line;
    @info=split(/\t/,$line);
    if($info[5]<0){
        if(exists $same{$info[7]}){
            $same{$info[7]}.=$info[2].'+'.$info[3].';';
            if($info[6]<$sp{$info[7]}){
                $sp{$info[7]}=$info[6];
            }
        }
        else{
            $same{$info[7]}=$info[2].'+'.$info[3].';';
            $sp{$info[7]}=$info[6];
        }
    }
}
=pod
$head{19}=1;$tail{19}=87;
$head{2}=88;$tail{2}=441;
$head{20}=442;$tail{20}=547;
$head{21}=548;$tail{21}=586;
$head{22}=587;$tail{22}=637;
$head{3}=638;$tail{3}=846;
$head{5}=847;$tail{5}=1090;
$head{6}=1091;$tail{6}=1343;
$head{1}=1344;$tail{1}=1711;
$head{7}=1712;$tail{7}=1918;
$head{8}=1919;$tail{8}=2126;
$head{9}=2127;$tail{9}=2295;
$head{23}=2296;$tail{23}=2384;
$head{10}=2385;$tail{10}=2577;
$head{11}=2578;$tail{11}=2748;
$head{12}=2749;$tail{12}=2968;
$head{13}=2969;$tail{13}=3084;
$head{14}=3085;$tail{14}=3211;
$head{15}=3212;$tail{15}=3350;
$head{16}=3351;$tail{16}=3435;
$head{17}=3436;$tail{17}=3561;
=cut
$time=1;
$m=0;
open REF,"./Knee_ref.txt";
while($line=<REF>){
    chomp $line;
    @info=split(/\t/,$line);
    if($m == 0){
        $rem=$info[0];
        $start=$time;
    }
    if($info[0] eq 'X'){$info0=23;} else{$info0=$info[0];}
    if($info0 ne $rem){
        $head{$rem}=$start;
        $tail{$rem}=$time-1;
        $start=$time;
        $rem=$info0;
    }
    $time++;
    $m=1;
    $tadh{$info[3]}=$info[1];
    $tadt{$info[3]}=$info[2];
    $chr{$info[3]}=$info[0];
}
$head{$rem}=$start;
$tail{$rem}=$time-1;
$line=<DIF>;
while($line=<DIF>){
    chomp $line;
    @info=split(/\t/,$line);
    if($info[5]>0){
        if((exists $same{$info[7]})||(exists $same{$info[8]})){
            if($info[5]>0){
                for($i=1;$i<=23;$i++){
                    if(($info[7]>=$head{$i} && $info[7]<=$tail{$i})&&($info[8]>=$head{$i} && $info[8]<=$tail{$i})){
                        $id1=$info[7].'_'.$info[8];
                        $id2=$info[8].'_'.$info[7];
                        if(exists $dif{$id1}){
                            $dif{$id1}.=$info[2].'+'.$info[3].';';
                            if($info[6]<$dp{$id1}){
                                $dp{$id1}=$info[6];
                            }
                        }
                        elsif(exists $dif{$id2}){
                            $dif{$id2}.=$info[2].'+'.$info[3].';';
                            if($info[6]<$dp{$id2}){
                                $dp{$id2}=$info[6];
                            }
                        }
                        elsif(!(exists $dif{$id1}) && !(exists $dif{$id2})){
                            $dif{$id1}=$info[2].'+'.$info[3].';';
                            $dp{$id1}=$info[6];
                        }
                    }
                }
            }
        }
    }
}
open DYS,">./corr_dys_bdry_flt.txt";
$j=1;
print DYS "no\tcorrelation_decreased_same_domain_pairs\tcorrelation_increased_cross_boundary_pairs\tbest_p_value\tdysregulated boundary\tTAD1\tTAD2\n";
foreach my $key(sort keys %dif){
    ($tad1,$tad2)=split(/_/,$key);
    if(exists $same{$tad1} && exists $same{$tad2}){
        $same=$same{$tad1}.$same{$tad2};
        if($sp{$tad1}>$sp{$tad2}){
            $sp=$sp{$tad2};
        }
        else{
            $sp=$sp{$tad1};
        }
    }
    elsif(exists $same{$tad1}){
        $same=$same{$tad1};
        $sp=$sp{$tad1};
    }
    elsif(exists $same{$tad2}){
        $same=$same{$tad2};
        $sp=$sp{$tad2};
    }
    else{
        next;
    }
    if($sp>$dp{$key}){
        $pv=$dp{$key};
    }
    else{
        $pv=$sp;
    }
    if(min($tadt{$tad1},$tadt{$tad2})-max($tadh{$tad1},$tadh{$tad2})<=0){
        if($tadh{$tad1}>=$tadh{$tad2}){
            $bdry=$chr{$tad1}.':'.$tadt{$tad2}.'-'.$tadh{$tad1};
        }
        elsif($tadh{$tad1}<$tadh{$tad2}){
            $bdry=$chr{$tad1}.':'.$tadt{$tad1}.'-'.$tadh{$tad2};
        }
        print DYS "$j\t$same\t$dif{$key}\t$pv\t$bdry\t$tadh{$tad1}:$tadt{$tad1}\t$tadh{$tad2}:$tadt{$tad2}\n";
        $j++;
    }
}
